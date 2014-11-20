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
with open("KEGG_Pathway_Genes.txt", "r") as fh:
    fh.readline().rstrip() # remove header
    
    kegg = [] # empty list to put KEGG Pathway Genes file in
    
    for line in fh.readlines(): # loop through all lines of file
        words = line.rstrip().split("\t") # split line by tab
        kegg.append(words) # add to master KEGG pathway list

