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

#########################
### GET PATHWAY GENES ###
#########################

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

# parse through HTML results
