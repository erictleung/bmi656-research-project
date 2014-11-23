#!/usr/bin/python

# Name:        Eric Leung
# Assignment:  Research Project
# Date:        December 8th, 2014
# script:      crossSpecies.py
# Description: Examine cross-species conservation of chosen pathway genes
# Example:     "python crossSpecies.py PATHWAY_NAME"

###################
### GET PATHWAY ###
###################

import sys
pathway = "\"" + sys.argv[1] + "\"" # put pathway string together

#############################
### EXTRACT PATHWAY GENES ###
#############################

# search for pathway ID
import urllib
opener = urllib.FancyURLopener({})
searchSite = "http://rest.kegg.jp/find/pathway/" # original site
straight = opener.open(searchSite+pathway) # put together query
rawResults = straight.read().rstrip().split("\t") # get results
pathID = rawResults[0].split(":")[1] # get pathway ID

# convert ID to human's
import re # import regex package
ID = re.findall(r"[0-9]+", pathID) # extract ID number
humanID = "hsa" + ID[0] # convert path ID to human

# get pathway genes
pathwaySite = "http://rest.kegg.jp/get/" # searching site
pathwayInfo = opener.open(pathwaySite + humanID)
rawResults = pathwayInfo.read() # get results
listOfResults = rawResults.split("\n") # put each line into list

genes = [] # empty liist to put all genes in

# parse through HTML page list
for line in listOfResults:
    temp = [] # list to put results
    words = line.strip().split() # remove white space
    if len(words) == 0: # if there is an empty line
        continue # skip the empty line
    elif words[0] == "GENE": # for the case with GENE first
        genes.append(words[2].strip(";")) # add gene name list
    elif re.match(r"\d+", words[0]): # ID number in front of gene
        entry = line.strip().split(";") # split into two parts

        # get ID numbere and name
        first = entry[0] # take first entry
        idName = first.split() # split string by whitespace
        temp.extend(idName) # add ID and gene name to temp

        # get gene full name
        second = entry[1] # take second entry
        name = second.split("[") # split by "[" char
        temp.append(name[0].strip())

        # put Gene ID, Name, Full name into final list
        genes.append(temp)

#############################
### OBTAIN GENE SEQUENCES ###
#############################

from Bio import Entrez

# provide email address
email = "leunge@ohsu.edu"
Entrez.email = email

handle = Entrez.esearch(db="protein", term=genes[0][2])
record = Entrez.read(handle)
ids = record["IdList"]
handle.close()
#for i in ids:
#    handle = Entrez.esummary(db="protein", id=i)
#    summary = Entrez.read(handle)
#    print summary
#    handle.close()

#########################
### CLUSTAL ALIGNMENT ###
#########################

from Bio.Align.Applications import ClustalwCommandline

# create Clustalw command
# command = ClustalwCommandline("clustalw2", infile = "")
