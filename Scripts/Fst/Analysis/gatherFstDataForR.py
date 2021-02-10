"""
This script is designed to collect Fst values, chromo, and coord for each hybrid 
group vs one of it's parent races and output the data in this R friendly format:
markerName FstHybCP  FstHybCB  FstHybSG chromo coord
Created by David E. Hufnagel
"""
import sys

fstCPf = open(sys.argv[1])     #fsts of HybCP vs one parent
fstCBf = open(sys.argv[2])     #fsts of HybCB vs one parent
fstSGf = open(sys.argv[3])     #fsts of HybSG vs one parent
markerData = open(sys.argv[4]) #Germplasm_survey-SNP_list_sorted.txt.mod
out = open(sys.argv[5], "w")   #parvVhybFsts.txt  OR  mexVhybFsts.txt


#Go through fstFD and make a dict of key: markerName val: Fst where negative fst values are set to 0
def GatherFstData(fstFD):
    fstDict = {}
    fstFD.readline()
    for line in fstFD:
        lineLst = line.strip().split()
        name = lineLst[0]
        fst = lineLst[1]
        
        #for fst=NA, do nothing
        if fst != "NA":
            fst = float(fst)
            if fst <0.0:
                fst = 0.0
                
            fstDict[name] = fst
            
    return(fstDict)
    
#Gather Fst data from all 3 files
fstDictCP = GatherFstData(fstCPf)
fstDictCB = GatherFstData(fstCBf)
fstDictSG = GatherFstData(fstSGf)

#Go through markerData and make a dict of key: markerName   val: (chromo, coordinate)
coordDict = {}
markerData.readline() #skip the first line
for line in markerData:
    lineLst = line.strip().split()
    name = lineLst[0]
    chromo = lineLst[2]
    coord = lineLst[3].strip('"').replace(",","")
    coordDict[name] = (chromo, coord)

#Output title line
titleLine = "markerName\tfstCP\tfstCB\tfstSG\tchromo\tcoord\n" % ()
out.write(titleLine)

#Go through all markers and output data from the 3 dictionaries
for name in fstDictCP.keys(): #CP was chosen arbitrarily
    printme = True
    fstCP = fstDictCP[name]
    
    if name in fstDictCB:
        fstCB = fstDictCB[name]
    else:
        printme = False
        
    if name in fstDictSG:
        fstSG = fstDictSG[name]
    else:
        printme = False
    
    if printme == True:
        chromo = coordDict[name][0]    
        coord = coordDict[name][1]
        
        if not chromo in ["unknown","multiple"] and not coord in ["unknown","multiple"]:
            newLine = "%s\t%s\t%s\t%s\t%s\t%s\n" % (name, fstCP, fstCB, fstSG, chromo, coord)
            out.write(newLine)


fstCPf.close()
fstCBf.close()
fstSGf.close()
markerData.close()
out.close()