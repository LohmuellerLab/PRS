import sys
import gzip

with gzip.open(sys.argv[1], 'rt') as f:
    for line in f:
        l = line.split("\t")
        if l[0] == "variant":
            continue
        else:
            pos_info = l[0].split(":")
            print("chr"+pos_info[0], int(pos_info[1])-1, pos_info[1], *l[1:], sep="\t", end="")
