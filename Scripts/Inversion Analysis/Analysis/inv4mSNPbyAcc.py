"""
This script is intended to take an input big info file, calculate PZD00030.1 
genotype frequencies, group by accesion, and create a file with the format:
accesion   hybStatus   altitude   PZD00030.1_cFreq   PZD00030.1_tFreq   PZD00030.1_nFreq
Created by David E. Hufnagel
"""
import sys
inp = open(sys.argv[1])      #ZeaAllInfo.pmz
out = open(sys.argv[2], "w") #inv4mSNPbyAcc.txt


def DetermineHybStatComplex(group):
    #Count each type per individual
    phcCnt = 0 #parvHC count
    mhcCnt = 0
    pamCnt = 0 #parv ambiguous count
    mamCnt = 0
    hybCnt = 0
    for pair in group:
        if pair[0] == "parvHC":
            phcCnt += 1 #parvHC
        elif pair[0] == "mexHC":
            mhcCnt += 1 #mexHC
        elif pair[1] == "Zea_mays_parviglumis" and pair[0] == "other":
            pamCnt += 1 #ambig parv
        elif pair[1] == "Zea_mays_mexicana" and pair[0] == "other":
            mamCnt += 1 #ambig mex
        elif pair[0] == "hybrid":
            hybCnt += 1 #hybrid
        else:
            print "ERROR"
            sys.exit()
    
    #Determine the majority group and the total minority percentage
    top = 0
    ind = 0
    bestInd = -9
    counts = [phcCnt, mhcCnt, pamCnt, mamCnt, hybCnt]
    for cnt in counts:
        if cnt > top:
            top = cnt
            bestInd = ind
        ind += 1
            
    majority = counts[bestInd]
    minority = sum(counts) - majority    
    minorFrac = float(minority) / float(sum(counts))

    #If minor faction .2 or higher, call the population mixed, otherwise call 
    #it by its majority
    if minorFrac >= 0.2:
        return "black"
    else:
        if bestInd == 0: #parvHC
            return "blue2"
        elif bestInd == 1: #mexHC
            return "firebrick4"
        elif bestInd == 2: #parv ambig
            return "cyan"
        elif bestInd == 3: #mex ambig
            return "firebrick1"
        elif bestInd == 4: #hybrids
            return "yellow1"
        else:
            print "ERROR1!"
            sys.exit() 

def SaveIntoDict(key,val,dictx):
    if key not in dictx:
        dictx[key] = [val,]
    else:
        dictx[key].append(val)
        
def CalcAlleleFreqs(genos):
    allAlleles = []
    for pair in genos:
        allAlleles.append(pair[0])
        allAlleles.append(pair[2])
    cfreq = float(allAlleles.count("C")) / float(len(allAlleles))
    tfreq = float(allAlleles.count("T")) / float(len(allAlleles))
    nfreq = float(allAlleles.count("N")) / float(len(allAlleles))
        
    return cfreq,tfreq, nfreq

#Go through title and determine the index of PZD00030.1
inp.readline()
markers = inp.readline().strip().split(":")[-1].split(",")
ind = markers.index("PZD00030.1")

#Go through inp and determine the genotype at PZD00030.1 for each ind and 
#  store in a dict of key: acc   val: genos as well as a
#  dict of key: acc val: [(hybA,taxA),(hybA,taxB),...]
accDict = {}
hybDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        acc = lineLst[3]
        genos = lineLst[-1].split(",")
        geno = genos[ind]
        hyb = lineLst[4]
        tax = lineLst[1]    
        
        if hyb != "NA":
            SaveIntoDict(acc,geno,accDict)

        if acc in hybDict:
            hybDict[acc].append((hyb,tax))
        else:
            hybDict[acc] = [(hyb, tax),]
            
#Go through accDict, calculate genetpye frequencies, and put them in accDict2 with key: acc val: (cfreq,tfreq,nfreq)
accDict2 = {}
for acc,genos in accDict.items():
    cfreq, tfreq, nfreq = CalcAlleleFreqs(genos)
    accDict2[acc] = (cfreq, tfreq, nfreq)

#Go through hybDict and make a determination of hybrid status for each acc
hybStatDict = {}
for acc, group in hybDict.items():
    #this means all the hybrid statuses are the same
    if len(set(group)) == 1:
        if group[0][0] == "parvHC":
            hybStat = "blue2" #parvHC
        elif group[0][0] == "mexHC":
            hybStat = "firebrick4" #mexHC
        elif group[0][1] == "Zea_mays_parviglumis" and group[0][0] == "other":
            hybStat = "cyan" #ambig parv
        elif group[0][1] == "Zea_mays_mexicana" and group[0][0] == "other":
            hybStat = "firebrick1" #ambig mex
        elif group[0][0] == "hybrid":
            hybStat = "darkorange2" #hybrid
        else:
            if group[0][0] != "NA":
                print "ERROR2"
                sys.exit()
            else:
                hybStat = "NA"
    elif len(group) == 2:
        hybStat = "black" #mixed
    else:
        hybStat = DetermineHybStatComplex(group)
    
    if hybStat != "NA":
        hybStatDict[acc] = hybStat

#Go through inp again and make a dict of key: acc, val: alt
inp.seek(0)
altDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        acc = lineLst[3]
        alt = lineLst[10]
        altDict[acc] = alt # I tested it, and all members of an accesion have the same altitude


#output results
newLine = "accesion\thybStatus\taltitude\tPZD00030.1_cFreq\tPZD00030.1_tFreq\tPZD00030.1_nFreq\n"
out.write(newLine)
for acc in accDict.keys():
    newLine = "%s\t%s\t%s\t%s\t%s\t%s\n" % (acc, hybStatDict[acc], altDict[acc],accDict2[acc][0],accDict2[acc][1],accDict2[acc][2])
    out.write(newLine)


inp.close()
out.close()
