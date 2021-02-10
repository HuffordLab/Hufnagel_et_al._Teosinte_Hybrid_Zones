In "Analysis", "makeHierFstat.py" prepares data for splitting by "splitHierFstat.py" and after locus-by-locus Fst calculation in "obtainingFST.R".

In "Analysis", "splitHierFstat.py" splits the one big hierFstat file created by "makeHierFstat.py" into one file for each pairwise comparison, which is ran through "obtainingFST.R" to calculate locus-by-locus Fst for all combinations of each hybrid group and each major race of parviglumis and mexicana.

In "Analysis", "obtainingFST.R" calculates locus-by-locus Fst for all combinations of each hybrid group and each major race of parviglumis and mexicana.

In "Analysis", "averageFsts.py" locus-by-locus Fst values are averaged to produce global Fst values.

In "Analysis", "gatherFstDataForR.py" collects locus-by-locus Fst values and chromosome and coordinate data for each marker for each hybrid group vs one of it's
 putative parent races and outputs the data in an R friendly format for use in "makeFstByFstFigs.R" and "makeFstByFstSubset.R".


In "Figure Generation", "makeFstByFstFigs.R" plots FstxFst figures for all chromosomes and both parviglumis and mexicana alleles. 

In "Figure Generation", "makeFstByFstSubset.R" plots FstxFst figures for select chromosomes and only mexicana alleles.

In "Figure Generation", "makeFstDistFrom45figs.R" plots distance from 45 degree figures based off of FstxFst figures for all chromosomes and both parviglumis and mexicana alleles.

In "Figure Generation", "makeFstDist45subset.R" plots distance from 45 degree
figures based off of FstxFst figures for select chromosomes and only mexicana alleles.
