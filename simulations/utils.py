import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def get_beta(s, tau):
    if tau == 0:
        C=0.10
    elif tau == 0.25:
        C=0.6
    elif tau == 0.5:
        C=1.8
    delta = np.random.choice([-1,1])
    eta = np.random.normal(0, 0.5)
    return(delta*np.power(abs(float(s)), tau)*C*(1+eta))

def score(genome_string, _betas, pop, id_map1, id_map2, freq1, freq2):
    mutations = [int(i) for i in genome_string.split(" ")[2:]]
    curr_betas = []
    for mut in mutations:
        if (pop == 1):
            _id = id_map1[mut]
        else:
            _id = id_map2[mut]
        try:
            curr_betas.append(_betas[_id])
        except KeyError:
            continue
    return(sum(curr_betas))

def score_shared(genome_string, _betas, pop, id_map1, id_map2, freq1, freq2):
    mutations = [int(i) for i in genome_string.split(" ")[2:]]
    curr_betas = []
    for mut in mutations:
        if (pop == 1):
            _id = id_map1[mut]
        else:
            _id = id_map2[mut]
        try:
            AF1 = freq1[_id]
        except KeyError:
            AF1 = 0
        try:
            AF2 = freq2[_id]
        except KeyError:
            AF2 = 0
        if AF1 > 0 and AF2 > 0:
            try:
                curr_betas.append(_betas[_id])
            except KeyError:
                curr_betas.append(0)
    return(sum(curr_betas))

def score_private(genome_string, _betas, pop, id_map1, id_map2, freq1, freq2):
    mutations = [int(i) for i in genome_string.split(" ")[2:]]
    curr_betas = []
    for mut in mutations:
        if (pop == 1):
            _id = id_map1[mut]
        else:
            _id = id_map2[mut]
        try:
            AF1 = freq1[_id]
        except KeyError:
            AF1 = 0
        try:
            AF2 = freq2[_id]
        except KeyError:
            AF2 = 0
        if pop == 2 and AF1 == 0 and AF2 > 0:
            try:
                curr_betas.append(_betas[_id])
            except KeyError:
                curr_betas.append(0)
        if pop == 1 and AF2 == 0 and AF1 > 0:
            try:
                curr_betas.append(_betas[_id])
            except KeyError:
                curr_betas.append(0)
    return(sum(curr_betas))


def parse(file, tau, prs=False):
    N_HAPLOID=10000
    betas1 = {}
    sel1 = {}
    freq1 = {}

    betas2 = {}
    sel2 = {}
    freq2 = {}

    id_map1 = {}
    id_map2 = {}

    pop1_info = {}
    pop2_info = {}

    prs_all1 = []
    prs_all2 = []

    prs_shared1 = []
    prs_shared2 = []
    prs_private1 = []
    prs_private2 = []

    afr_genomes = []

    muts = 0
    genomes = 0
    with open(file, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if "#OUT:" in line:
                continue
            if "Mutations" in line:
                muts += 1
                genomes += 1
                continue
            if "Genomes" in line:
                genomes += 1
                muts += 1
                continue
            if muts == 1:
                id_map1[int(line.split(" ")[0])] = int(line.split(" ")[1])
                sel1[int(line.split(" ")[1])] = float(line.split(" ")[4])
                betas1[int(line.split(" ")[1])] = get_beta(line.split(" ")[4], tau)
                freq1[int(line.split(" ")[1])] = int(line.split(" ")[8])/N_HAPLOID
                pop1_info[int(line.split(" ")[1])] = [float(line.split(" ")[4]), get_beta(line.split(" ")[4], tau), int(line.split(" ")[8])/N_HAPLOID]
            if muts == 3:
                id_map2[int(line.split(" ")[0])] = int(line.split(" ")[1])
                sel2[int(line.split(" ")[1])] = float(line.split(" ")[4])
                betas2[int(line.split(" ")[1])] = get_beta(line.split(" ")[4], tau)
                freq2[int(line.split(" ")[1])] =  int(line.split(" ")[8])/N_HAPLOID
                pop2_info[int(line.split(" ")[1])] = [float(line.split(" ")[4]), get_beta(line.split(" ")[4], tau), int(line.split(" ")[8])/N_HAPLOID]

            if genomes == 2:
                afr_genomes.append(line)
            if genomes == 4:
                prs_all2.append(score(line, betas2, 2, id_map1, id_map2, freq1, freq2))
                prs_shared2.append(score_shared(line, betas2, 2, id_map1, id_map2, freq1, freq2))
                prs_private2.append(score_private(line, betas2, 2, id_map1, id_map2, freq1, freq2))

        for g in afr_genomes:
            prs_all1.append(score(g, betas1, 1, id_map1, id_map2, freq1, freq2))
            prs_shared1.append(score_shared(g, betas1, 1, id_map1, id_map2, freq1, freq2))
            prs_private1.append(score_private(g, betas1, 1, id_map1, id_map2, freq1, freq2))

        p1 = pd.DataFrame.from_dict(pop1_info, orient="index")
        p2 = pd.DataFrame.from_dict(pop2_info, orient="index")
        merged = p1.merge(p2, left_index=True, right_index=True, how='outer')
        merged.columns = ["p1.s", 'p1.b', 'p1.f', 'p2.s', 'p2.b', 'p2.f']
        prs_dict = {
            'all1':prs_all1, "private1":prs_private1, "shared1":prs_shared1,
            'all2':prs_all2, "private2":prs_private2, "shared2":prs_shared2}
        prs_df = pd.DataFrame.from_dict(prs_dict)
        if prs:
            return(prs_df)
        else:
            return(merged)
