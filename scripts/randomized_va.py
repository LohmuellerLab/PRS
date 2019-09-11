import gzip
from itertools import islice
import numpy as np
import math
import sys
import argparse

parser = argparse.ArgumentParser(description="Script to calculate Va from summary stat data")
parser.add_argument("-w", "--window_size", action="store", required=True, help="Window size in terms of number of SNPs")
parser.add_argument("-m", "--sample_size", action="store", required=True, help="Number of SNPs to resample per window")
parser.add_argument("-f", "--file", action="store", required=True, help="GZIP'd input file")
parser.add_argument("-a", "--af", action="store", required=True, help="Allele frequency cutoff for calling SNPs private")
parser.add_argument("-t", "--thresh", action="store", default=1, required=False,help="P value cutoff for SNPs")
args = parser.parse_args()


W = int(args.window_size)#50000 # window size in terms of number of SNPs
M = int(args.sample_size) # number of snps to sample in a window
header=True

def jk_resamp(data):
    n = data.shape[0]
    if n <= 0:
        raise ValueError("data must contain at least one measurement.")

    resamples = np.empty([data.shape[0], data.shape[0]-1, data.shape[1]])

    for i in range(n):
        resamples[i] = np.delete(data, i, 0)

    return resamples

def jk_stat(data, statistic, conf_lvl=0.95):
    from scipy.special import erfinv

    # make sure original data is proper
    n = data.shape[0]
    if n <= 0:
        raise ValueError("data must contain at least one measurement.")

    resamples = jk_resamp(data)

    stat_data = statistic(data)
    jack_stat = []
    for i in resamples:
        jack_stat.append(np.sum(i[:,0])/sum(i[:,1]))
    jack_stat = np.array(jack_stat)
    mean_jack_stat = np.mean(jack_stat, axis=0)

    # jackknife bias
    bias = (n-1)*(mean_jack_stat - stat_data)

    # jackknife standard error
    std_err = np.sqrt((n-1)*np.mean((jack_stat - mean_jack_stat)*(jack_stat -
                                    mean_jack_stat), axis=0))

    # bias-corrected "jackknifed estimate"
    estimate = stat_data - bias

    # jackknife confidence interval
    if not (0 < conf_lvl < 1):
        raise ValueError("confidence level must be in (0, 1).")

    z_score = np.sqrt(2.0)*erfinv(conf_lvl)
    conf_interval = estimate + z_score*np.array((-std_err, std_err))

    return estimate, bias, std_err, conf_interval




big_list = []
with gzip.open(args.file,'rt') as infile:
    if(header):
        next(infile)
    while True:
        lines_gen = islice(infile, W)
        if not lines_gen:
            break
        try:
            l = np.random.choice(list(lines_gen), M)
        except ValueError:
            break
        va_w = []
        va_pw = []
        for snp in l:
            s = snp.split("\t")
            AF = float(s[2])
            pval = float(s[-1])
            if (pval > float(args.thresh)):
                continue
            if len(s) == 12:
                beta = float(s[8])
            elif len(s) == 11:
                beta = float(s[7])
            if math.isnan(AF) or math.isnan(beta):
                continue
            va_w.append(2*AF*(1-AF)*beta**2)
            if (AF <= float(args.af)):
                va_pw.append(2*AF*(1-AF)*beta**2)

        try:
            va_p = sum(va_pw)
            va = sum(va_w)
        except ZeroDivisionError:
            continue
        if va_p > va:
            va = va_p # constrain va !> va_p
        big_list.append([va_p,va])

test_statistic = lambda x: np.sum(x[:,0])/sum(x[:,1])
a = np.array(big_list)
estimate, bias, stderr, conf_interval = jk_stat(a, test_statistic)
print(estimate, conf_interval[0], conf_interval[1], stderr, bias, W, M, args.af, args.thresh, args.file)
