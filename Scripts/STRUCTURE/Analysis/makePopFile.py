#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is designed to take the big info file which contains individual
data and use it to create a population file which contains population data
with the format:
popID   taxonomy   hybrid_status   hybrid_group   latitude   longitude   country   state   altitude   pop_avg_STRUCTURE_q-matrix_vals(Parv_Mex_Maize)
Created by David E. Hufnagel on Tue Dec 22, 2020
"""
import sys
inp = open(sys.argv[1])
out = open(sys.argv[2], "w")



def SaveIntoDict(key, val, dictx):
    if not key in dictx:
        dictx[key] = [val,]
    else:
        dictx[key].append(val)
        
def Avg(listx):
    avg = float(sum(listx)) / float(len(listx))
    return(avg)



#Go through inp and make a dict of key: popID val: (tax,hyb,hybGroup,lat,lon,cntry,state,alt,qs) for hybrids
popDict = {}
inp.readline(); inp.readline()
for line in inp:
    lineLst = line.strip().split("\t")
    tax = lineLst[1]; hyb = lineLst[4]; hybGroup = lineLst[5]; lat = lineLst[6]
    lon = lineLst[7]; cntry = lineLst[8]; state = lineLst[9]; alt = lineLst[10]
    qs = lineLst[11]
    
    if hyb == "hybrid":
        popID = "%s_%s" % (lat,lon)
        val = (tax,hyb,hybGroup,lat,lon,cntry,state,alt,qs)
        SaveIntoDict(popID, val, popDict)
        
        
#Create the output title
##Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))
out.write("#popID\ttaxonomy\thybrid(parvHC,mexHC,hybrid,other)\thybrid_group(Central_Plateau,Central_Balsas,South_Guerrero)\tlatitude\tlongitude\tcountry\tstate\taltitude\tpop_avg_STRUCTURE_q-matrix_vals(Parv_Mex_Maize)\n")
        
    
#Go through dict, combine information and output the result
popNum = 1
for popID, data in popDict.items():
    if len(data) == 1:
        tax = data[0][0]; hyb = data[0][1]; hybGroup = data[0][2]
        lat = data[0][3]; lon = data[0][4]; cntry = data[0][5]
        state = data[0][6]; alt = data[0][7]; qs = data[0][8]
        newline = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
            (popNum, tax, hyb, hybGroup, lat, lon, cntry, state, alt, qs)
    else:
        taxs = []; hybs = []; hybGroups = []; lats = []; lons = []; cntrys = []
        states = []; alts = []; qss = []
        for i in range(len(data)):
            taxs.append(data[i][0]); hybs.append(data[i][1])
            hybGroups.append(data[i][2]); lats.append(data[i][3])
            lons.append(data[i][4]); cntrys.append(data[i][5])
            states.append(data[i][6]); alts.append(data[i][7])
            qss.append(data[i][8])
            
        if len(set(taxs)) != 1:
            print(taxs)
            sys.exit()
        else:
            tax = taxs[0]
            
        if len(set(hybs)) != 1:
            print(hybs)
            sys.exit()
        else:
            hyb = hybs[0]
            
        if len(set(hybGroups)) != 1:
            print(hybGroups)
            sys.exit()
        else:
            hybGroup = hybGroups[0]
            
        if len(set(lats)) != 1:
            print(lats)
            sys.exit()
        else:
            lat = lats[0]
            
        if len(set(lons)) != 1:
            print(lons)
            sys.exit()
        else:
            lon = lons[0]
            
        if len(set(cntrys)) != 1:
            print(cntrys)
            sys.exit()
        else:
            cntry = cntrys[0]
            
        if len(set(states)) != 1:
            print(states)
            sys.exit()
        else:
            state = states[0]
            
        if len(set(alts)) != 1:
            alt = "NA"
        else:
            alt = alts[0]
            
        if len(set(qss)) != 1:
            ps = []; ms = []; zs = []  #These stand for parvs, mexs, and maizes
            for temp in qss:
                tempLst = temp.split("_")
                ps.append(float(tempLst[0]))
                ms.append(float(tempLst[1]))
                zs.append(float(tempLst[2]))
            avgP = Avg(ps); avgM = Avg(ms); avgZ = Avg(zs)
            qs = "%s_%s_%s" % (avgP, avgM, avgZ)
        else:
            qs = qss[0]

    
        newline = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
            (popNum, tax, hyb, hybGroup, lat, lon, cntry, state, alt, qs)
    
    out.write(newline)
    popNum += 1




inp.close()
out.close()