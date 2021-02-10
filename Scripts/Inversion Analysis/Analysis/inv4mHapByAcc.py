"""
This script is designed to take as input a big data file, calculate inv4m 
haplotype counts, group by accesion, and create an
output of: accesion   hybStatus   altitude   inv4mPcount   inv4mMcount   inv4mHcount   inv4mNcount   inv4mRcount
Created by David E. Hufnagel
"""
import sys
inp = open(sys.argv[1])      #ZeaAllInfo.pmz
out = open(sys.argv[2], "w") #inv4mHapByAcc.txt


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
            print "ERROR1"
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
            print "ERROR2!"
            sys.exit()   
     
    
def SaveIntoDict(key,val,dictx):
    if key not in dictx:
        dictx[key] = [0,0,0,0,0]
        
    if val == "P":
        dictx[key][0] += 1
    elif val == "M":
        dictx[key][1] += 1
    elif val == "H":
        dictx[key][2] += 1
    elif val == "N":
        dictx[key][3] += 1
    elif val == "R":
        dictx[key][4] += 1
    else:
        print "ERROR4!"
        sys.exit()
        
#Go thrugh inp and make a dict of key: acc and  val: inv4m status as well as a
#  dict of key: acc val: [(hybA,taxA),(hybA,taxB),...]
accDict = {}
hybDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        acc = lineLst[3]; inv = lineLst[12]; hyb = lineLst[4]; tax = lineLst[1]
        
        if hyb != "NA":
            SaveIntoDict(acc,inv,accDict) 
            
            if acc in hybDict:
                hybDict[acc].append((hyb,tax))
            else:
                hybDict[acc] = [(hyb, tax),]

#Go through hybDict and make a determination of hybrid status for each group
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
            print "ERROR3"
            sys.exit()
    elif len(group) == 2:
        hybStat = "black" #mixed
    else:
        hybStat = DetermineHybStatComplex(group)
    
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
newLine = "accesion\thybStatus\taltitude\tinv4mPcount\tinv4mMcount\tinv4mHcount\tinv4mNcount\tinv4mRcount\n"
out.write(newLine)
for acc in accDict.keys():
    newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (acc, hybStatDict[acc], altDict[acc], accDict[acc][0],accDict[acc][1],accDict[acc][2],accDict[acc][3],accDict[acc][4])
    out.write(newLine)


inp.close()
out.close()
