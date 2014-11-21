#!/usr/bin/python

# Name: Eric Leung
# Assignment: Final Project
# Script: odds.py
# Description:
#    Use Python to perform odds ratio calculation to identify pathways
#    containing a larger number of differentially expressed genes than would
#    be expected by chance
# Example: "python odds.py"

##################
### FILE CHECK ###
##################

#  make sure the right files are in the current directory
import os.path, sys # import necessary packages
files = ["H5N1_VN1203_DE_Probes.txt",
         "H5N1_VN1203_UNIVERSE_Probes.txt",
         "KEGG_Pathway_Genes.txt"]
for n in files: # loop through list of files that should be directory
    temp = os.path.isfile(n) # check if file is in directory
    if temp != True: # if you get a false output from line above
        sys.exit("One of the following files for this script are not" + \
                 " in your current directory:" + "\n" + \
                 str(files)) # list of files

###################
### IMPORT DATA ###
###################

# focus on KEGG Pathway Genes
# SAMPLE FORMAT
# ID	PATHWAY TITLE	MEMBERS
# 04142	Lysosome	AP3S2	CTSE	GALNS   ...
# 04916	Melanogenesis	WNT10A	ADCY7	CREB3L3 ...
with open("KEGG_Pathway_Genes.txt", "r") as fh:
    fh.readline().rstrip() # remove header
    
    kegg = [] # empty list to put KEGG Pathway Genes file in
    
    for line in fh.readlines(): # loop through all lines of file
        words = line.rstrip().split("\t") # split line by tab
        kegg.append(words) # add to master KEGG pathway list

# focus on H5N1 VN1203 UNIVERSE Probes
# SAMPLE FORMAT
# A_23_P100011	AP3S2
# A_23_P100022	SV2B
# A_23_P100111	CHP
with open("H5N1_VN1203_UNIVERSE_Probes.txt", "r") as fh:
    fh.readline().rstrip() # remove header

    universe = {} # empty dictionary for universe probes

    for line in fh.readlines():
        words = line.rstrip().split("\t") # split line by tab
        universe[words[0]] = words[1] # save in dictionary

    universeProbes = universe.values() # get all probes in universe
    
# focus on H5N1 VN1203 DE Probes
# SAMPLE FORMAT
# A_23_P100539	ABCC6
# A_23_P100642	PNMT
# A_23_P100704	MAPK7
with open("H5N1_VN1203_DE_Probes.txt", "r") as fh:
    fh.readline().rstrip() # remove header

    de = {} # empty dictionary for universe probes

    for line in fh.readlines():
        words = line.rstrip().split("\t") # split line by tab
        de[words[0]] = words[1] # save in dictionary

    deProbes = de.values() # get all DE probes

#################
### ODD RATIO ###
#################

# variables/data types
sigPathways = [] # list of pathways with signif num DE genes

# matrix for odds ratio
#                | DE Genes | Non-DE 
# target pathway |    a     |    b 
# non-pathway    |    c     |    d

# loop through pathways to get genes
for pathway in kegg: # loop through all pathways
    name = pathway[0:2] # get pathway ID and pathway name
    genes = pathway[2:] # get probes in certain pathway
    genes = set(genes) # make a set of genes

    # make other data into sets
    universeSet = set(universeProbes)
    deSet = set(deProbes)

    # set values based on matrix shown above
    nonDE = universeSet.difference(deSet) # get non-DE genes
    a = float(len(genes.intersection(deSet)))
    b = float(len(genes.intersection(nonDE)))
    c = float(len(deSet.difference(genes)))
    d = float(len(nonDE.intersection(universeSet.difference(genes))))

    # find pathways with odds ratio greater than 1.5
    if (a*d)/(b*c) > 1.5:
        sigPathways.append(name)
