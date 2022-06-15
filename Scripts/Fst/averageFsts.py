"""
This script is designed to average local Fst values to produce a global Fst
value.
Created by David E. Hufnagel
"""
import sys

inp = open(sys.argv[1])

inp.readline()
fsts = []
cnt = 0 #a counter to keep track of non-NA values so they aren't included in the divisor
for line in inp:
    lineLst = line.strip().split()
    if lineLst[1] != "NA":
        fst = float(lineLst[1])
        
        #set negative fst values to 0
        if fst < 0:
            fst = 0
        fsts.append(fst)
        cnt += 1

print sum(fsts) / float(cnt)


inp.close()