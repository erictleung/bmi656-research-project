#!/usr/bin/python

# Name:        Eric Leung
# Assignment:  Research Project
# Date:        December 8th, 2014
# Python Ver:  2.7.8
# script:      crossSpecies.py
# Description: Examine cross-species conservation of chosen pathway genes
# Example:     "python crossSpecies.py PATHWAY_NAME"

#################
### FUNCTIONS ###
#################

def parse_kegg_html(pathway, species):
    """
    INPUT: pathway name, species of interest
    OUTPUT: HTML output from KEGG API
    """
    # search for pathway ID
    opener = urllib.FancyURLopener({})
    searchSite = "http://rest.kegg.jp/find/pathway/" # original site
    straight = opener.open(searchSite+pathway) # put together query
    rawResults = straight.read().rstrip().split("\t") # get results
    pathID = rawResults[0].split(":")[1] # get pathway ID

    # convert ID to human's
    ID = re.findall(r"[0-9]+", pathID) # extract ID number
    ID = species + ID[0] # convert path ID to human
    
    # get pathway genes
    pathwaySite = "http://rest.kegg.jp/get/" # searching site
    pathwayInfo = opener.open(pathwaySite + ID)
    rawResults = pathwayInfo.read() # get results
    listOfResults = rawResults.split("\n") # put each line into list
    
    return listOfResults

def get_genes(listOfResults):
    """
    INPUT: HTML output from KEGG API
    OUTPUT: gene ID, gene name, gene description in a list
    NOTE: need to run get_pathway function prior
    """
    genes = [] # empty liist to put all genes in
    
    # parse through HTML page list
    for line in listOfResults:
        temp = [] # list to put results
        words = line.strip().split() # remove white space
        if len(words) == 0: # if there is an empty line
            continue # skip the empty line
        elif words[0] == "GENE": # for the case with GENE first
            entry = line.strip().split(";") # split into two parts
    
            # get ID numbers
            first = entry[0] # take first entry
            idName = first.split() # split string by whitespace
            temp.extend(idName[1:]) # add ID and gene name to temp
    
            # get gene full name
            second = entry[1] # take second entry
            name = second.split("[") # split by "[" char
            temp.append(name[0].strip())
    
            # put Gene ID, Name, Full name into final list
            genes.append(temp)
        elif re.match(r"\d+", words[0]): # ID number in front of gene
            entry = line.strip().split(";") # split into two parts
    
            # get ID number and name
            first = entry[0] # take first entry
            idName = first.split() # split string by whitespace
            temp.extend(idName) # add ID and gene name to temp
    
            # get gene full name
            second = entry[1] # take second entry
            name = second.split("[") # split by "[" char
            temp.append(name[0].strip())
    
            # put Gene ID, Name, Full name into final list
            genes.append(temp)

    # return list of genes
    return genes

def get_accession(geneSearch):
    """
    INPUT: gene ID number
    OUTPUT: gene accession number
    """
    # create search
    handle = Entrez.efetch(db="gene", id=geneSearch, retmode="xml")
    fasta_record = handle.read() # XML file of gene
    handle.close() # close connection to database

    # parse through XML file
    root = et.fromstring(fasta_record)
    locus = root[0].find("Entrezgene_locus")
    product = locus[0].find("Gene-commentary_products")
    for access in product:
        name = access.find("Gene-commentary_accession").text
        if "NM" in name: # if it is a refseq nucleotide
            return name

###################
### GET PATHWAY ###
###################
"""
Take in argument from command line as target pathway
"""

import sys
pathway = "\"" + sys.argv[1] + "\"" # put pathway string together

#############################
### EXTRACT PATHWAY GENES ###
#############################
"""
Extract gene list for each species
"""

import re
import urllib

# species legend
species = {"Human":"hsa",
           "Mouse":"mmu",
           "Chimp":"ptr"}

# dictionary for each species
genes = {}

# loop through species
for org in species.keys():
    orgList = parse_kegg_html(pathway, species[org])
    genes[org] = get_genes(orgList)

#############################
### OBTAIN mRNA ACCESSION ###
#############################
"""
Get mRNA accession numbers for each gene
"""

from Bio import Entrez
import xml.etree.ElementTree as et

# provide email address
email = "leunge@ohsu.edu"
Entrez.email = email # add email to object

# dictionary of dictionaries for accession IDs for each species
allAccession = {}

# loop through all species
for org in genes.keys():
    # loop through genes
    for gene in genes[org]:
        target = gene[0] # focus on ID
        geneName = gene[1].upper() # focus on uppercase name
        if geneName not in allAccession.keys(): # if first instance
            allAccession[geneName] = {} # make dictionary
            allAccession[geneName][org] = get_accession(target)        

#############################
### OBTAIN mRNA SEQUENCES ###
#############################
"""
Loop through accession numbers to get sequences
"""

# dictionary of dictionaries for sequences for each species
allSequences = {}

#########################
### CLUSTAL ALIGNMENT ###
#########################

from Bio.Align.Applications import ClustalwCommandline

# create Clustalw command
# command = ClustalwCommandline("clustalw2", infile = "")
