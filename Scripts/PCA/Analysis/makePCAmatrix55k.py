"""
This script is designed to use the version 2 55k datafile to make a file for 
the purpose of making a PC plot in R.
Created By David E. Hufnagel on Sat May  5, 2018"""
import sys
inp = open(sys.argv[1])      #55kTeoNewSNPs_v4.txt
out = open(sys.argv[2], "w") #55kTeoNewSNPs_v4.txt.PCAmatrix.noMz

#Output header info
out.write("#python %s\n" % (" ".join(sys.argv)))

inp.readline(); inp.readline()


for line in inp:
    lineLst = line.strip().split("\t")   
    name = lineLst[0]
    out.write("%s\t" % (name))
    SNPs = lineLst[-1].split(",")
    cnt = 1
    for pair in SNPs:
        if cnt == len(SNPs):
            if pair == "NA":
                out.write("?/?\n")
            else:
                pair = "%s/%s\n" % (pair[0], pair[1])
                out.write(pair)
        else:
            if pair == "NA":
                out.write("?/?\t")
            else:
                pair = "%s/%s\t" % (pair[0], pair[1])
                out.write(pair)
        cnt += 1


inp.close()
out.close()
