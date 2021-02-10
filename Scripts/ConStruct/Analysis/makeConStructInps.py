"""
This script is designed to create 2/3 of the inputs for the R package conStruct:
the allele frequencies matrix and the coordinates matrix.  Individuals are
subsampled randomly at one individual per accesion to reduce sample size, and 
Maize is excluded. Individuals are ordered by group in this order: parvHC,
parvAmb,hybSG,hybCB,hybParvOther,hybCP,hybMexOther,mexAmb,mexHC
Created by David E. Hufnagel
"""
import sys, random

inp = open(sys.argv[1])               #ZeaAllInfo.pmz
alleleFreqFD = open(sys.argv[2], "w") #ZeaAllInfo.pmz.conS.alFreq
coordFD = open(sys.argv[3], "w")      #ZeaAllInfo.pmz.conS.coords


def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val) 
        

#Output title lines for both output files
inp.readline()
markers = inp.readline().strip().split(":")[-1].split(",")

alFreqTitle = "\t\t%s\n" % ("\t".join(markers))
alleleFreqFD.write(alFreqTitle)
coordTitle = "Lon\tLat\n"
coordFD.write(coordTitle)


#Go through inp, store genotype info in a dict of key: markerInd  val: all 
#  genotpyes at marker.  Also make a dict of key: accesion val: individuals 
#  for subsampling.  Also make a list of (group, name)
genoDict = {}
accDict = {}
groupLst = []
for line in inp:
    lineLst = line.strip().split("\t")
    name = lineLst[0]
    tax = lineLst[1]
    acc = lineLst[3]
    hyb = lineLst[4]
    pop = lineLst[5]
    genos = lineLst[-1].split(",")
    
    #Save into dictionary
    ind = 0
    for geno in genos:
        SaveIntoDict(ind, geno, genoDict) #These will stay ordered even in a dict because of the chosen key
        ind += 1
        
    #determine group for groupDict
    if tax == "Zea_mays_parviglumis":
        if hyb == "parvHC":
            group = "A_parvHC"
        elif hyb == "hybrid":
            if pop == "Central_Balsas":
                group = "D_hybCB"
            elif pop == "South_Guerrero":
                group = "C_hybSG"
            elif pop == "NA":
                group = "E_hybOther"
        elif hyb == "other":
            group = "B_parvAmb"
        else:
            print "error 1"
            sys.exit()
    elif tax == "Zea_mays_mexicana":
        if hyb == "mexHC":
            group = "I_mexHC"
        elif hyb == "hybrid":
            if pop == "Central_Plateau":
                group = "F_hybCP"
            elif pop == "NA":
                group = "G_hybOther"
        elif hyb == "other":
            group = "H_mexAmb"
        else:
            print "error 2"
            sys.exit()
    else:
        group = "NA"
    
    #build groupLst
    if tax != "Zea_mays_mays":
        SaveIntoDict(acc, name, accDict)
        groupLst.append((group, name))
    
#Go through groupLst and make a list of names in group order
groupLst.sort()
orderedNames = []
for pair in groupLst:
    orderedNames.append(pair[1])
     
#Subsample accDict randomly such that each acc has one individual
accDictSub = {} #accDict subsampled
for acc, inds in accDict.items():
    if len(inds) > 1:
        accDictSub[acc] = random.choice(inds)
    else:
        accDictSub[acc] = inds[0]
        
#Go through genoDict and determine the "1 state" arbitrarily, where idividuals 
#  who have a given genotype of homozygous at the "1 state" get a score of 1, 
#  homozygous alternative allele get's a score of 0, and heterozygotes get a 
#  score of 0.5.
oneStates = [0]*len(genoDict.keys())
for markerInd in genoDict.keys():
    for geno in genoDict[markerInd]:
        if geno[0] != "N":
            oneStates[markerInd] = geno[0]
            break #once a non-missing-data allele is found stop the loop

#Go through inp again, convert genotypes to "allele frequencies" for each 
#  individual based on "1 states". Also subsample 1 individual per accesion and
#  store information for both output files in dictionaries of key: name val: newLine
inp.seek(0)
alFreqNewLines = {} 
coordNewLines = {}  
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        name = lineLst[0]
        lat = lineLst[6]; lon = lineLst[7]
        genos = lineLst[-1].split(",")
        ind = 0
        alFreqs = []
        for geno in genos:
            oneState = oneStates[ind]
            alFreq = float(geno.count(oneState)) / 2.0
            alFreqs.append(str(alFreq))
            ind += 1
        
        if name in accDictSub.values():
            #output to alleleFreq file
            newLine = "%s\t%s\n" % (name, "\t".join(alFreqs))
            alFreqNewLines[name] = newLine
            
            #output to coordinate file
            newLine = "%s\t%s\n" % (lon,lat)
            coordNewLines[name] = newLine
            
#Use orderedNames and newLine dictionaries to print individuals in the desired 
#  order for both ouptuts
for name in orderedNames:
    if name in accDictSub.values():
        alleleFreqFD.write(alFreqNewLines[name])
        coordFD.write(coordNewLines[name])
    

inp.close()
alleleFreqFD.close()
coordFD.close()