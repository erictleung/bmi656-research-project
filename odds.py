#!/usr/bin/python

# Name:       Eric Leung
# Assignment: Final Project
# Date:       December 8th, 2014
# Python ver: 2.7.8
# Script:     odds.py
# Description:
#    Use Python to perform odds ratio calculation to identify pathways
#    containing a larger number of differentially expressed genes than would
#    be expected by chance
# Example:    "python odds.py"

##################
### FILE CHECK ###
##################

#  make sure the right files are in the current directory
import os.path, sys # import necessary packages
files = ["H5N1_VN1203_DE_Probes.txt",
         "H5N1_VN1203_UNIVERSE_Probes.txt",
         "KEGG_Pathway_Genes.txt"]
print "Checking for necessary files for analysis..."
for n in files: # loop through list of files that should be directory
    temp = os.path.isfile(n) # check if file is in directory
    if temp != True: # if you get a false output from line above
        sys.exit("One of the following files for this script are not" + \
                 " in your current directory:" + "\n" + \
                 str(files)) # list of files
print "Files necessary for analysis are in the current directory.\n"

###################
### IMPORT DATA ###
###################

# focus on KEGG Pathway Genes
# SAMPLE FORMAT
# ID	PATHWAY TITLE	MEMBERS
# 04142	Lysosome	AP3S2	CTSE	GALNS   ...
# 04916	Melanogenesis	WNT10A	ADCY7	CREB3L3 ...
print "Loading KEGG_Pathway_Genes.txt file..."
with open("KEGG_Pathway_Genes.txt", "r") as fh:
    fh.readline().rstrip() # remove header
    
    kegg = [] # empty list to put KEGG Pathway Genes file in
    
    for line in fh.readlines(): # loop through all lines of file
        words = line.rstrip().split("\t") # split line by tab
        kegg.append(words) # add to master KEGG pathway list
print "Successfully loaded KEGG_Pathway_Genes.txt\n"

# focus on H5N1 VN1203 UNIVERSE Probes
# SAMPLE FORMAT
# A_23_P100011	AP3S2
# A_23_P100022	SV2B
# A_23_P100111	CHP
print "Loading Universe probes..."
with open("H5N1_VN1203_UNIVERSE_Probes.txt", "r") as fh:
    fh.readline().rstrip() # remove header

    universe = {} # empty dictionary for universe probes

    for line in fh.readlines():
        words = line.rstrip().split("\t") # split line by tab
        universe[words[0]] = words[1] # save in dictionary

    universeGenes = universe.values() # get all genes in universe
    universeSet = set(universeGenes)
print "Successfully loaded Universe probes\n"
    
# focus on H5N1 VN1203 DE Probes
# SAMPLE FORMAT
# A_23_P100539	ABCC6
# A_23_P100642	PNMT
# A_23_P100704	MAPK7
print "Loading differentially expressed file..."
with open("H5N1_VN1203_DE_Probes.txt", "r") as fh:
    fh.readline().rstrip() # remove header

    de = {} # empty dictionary for universe probes

    for line in fh.readlines():
        words = line.rstrip().split("\t") # split line by tab
        de[words[0]] = words[1] # save in dictionary

    deGenes = de.values() # get all DE genes
    deSet = set(deGenes)
print "Successfully loaded differentially expressed data\n"

#################
### ODDS RATIO ###
#################

print "Calculating odds ratio for pathways...\n"

# variables/data types
sigPathways = [] # list of pathways with signif num DE genes

# matrix for odds ratio
#                | DE Genes | Non-DE 
# target pathway |    a     |    b 
# non-pathway    |    c     |    d

de_pathways = {} # dictionary to store DE and non-DE genes for pathway

# loop through pathways to get genes for target pathway
for pathway in kegg: # loop through all pathways
    name = pathway[0:2] # get pathway ID and pathway name
    print "Calculating odds ratio for " + name[1] # print pathway name
    genes = pathway[2:] # get probes in certain pathway
    print name[1] + " has " + str(len(genes)) + " number of genes."
    genes = set(genes) # make a set of genes in pathway

    # set values based on matrix shown above
    keggDE = universeSet.intersection(deSet) # make sure DE genes in KEGG
    nonDE = universeSet.difference(keggDE) # get non-DE genes

    # odds ratio pieces as defined above
    aSet = genes.intersection(keggDE)
    a = float(len(aSet))
    bSet = genes.intersection(nonDE)
    b = float(len(bSet))
    cSet = len(keggDE.difference(genes))
    c = float(cSet)
    dSet = len(nonDE.difference(genes))
    d = float(dSet)

    # add genes to pathway dictionary for DE and non-DE genes
    allList = [aSet, bSet] # [DE genes, non-DE genes]
    de_pathways[name[1]] = allList # put into dictionary

    # find pathways with odds ratio greater than 1.5
    OR  = (a*d)/(b*c)
    if OR > 1.5:
        sigPathways.append(name)
    print ""

####################
### SIG PATHWAYS ###
####################

print "Here are pathways with an OR >1.5:"
for path in sigPathways:
    print path[1]

######################
### CHOOSE PATHWAY ###
######################

# write to a file
fh = open("pathway_info.csv", "w")
myPathway = sigPathways[1][1] # select pathway to study
fh.write(myPathway + "\n") # write pathway name first
fh.write(",".join(de_pathways[myPathway][0]) + "\n") # DE genes
fh.write(",".join(de_pathways[myPathway][1]) + "\n") # non-DE genes
fh.close() # close file
