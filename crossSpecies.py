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
    temp = [] # line of accession numbers
    
    # get all mRNA accession numbers for the gene
    for access in product:
        name = access.find("Gene-commentary_accession").text
        #if "NM" in name: # if it is a refseq nucleotide
        #    return name
        temp.append(name)
    
    # return RefSeq mRNA if it exists
    for num in temp:
        if "NM" in num:
            print "Got " + num + "\n"
            return num

    # return predicted mRNA if it exists
    for num in temp:
        if "XM" in num:
            print "Got " + num + "\n"
            return num

    # return something
    print "Got " + num + "\n"
    return temp[0]

def keep_genes_common_with_humans(genes):
    """
    INPUT: genes list from each species
    OUTPUT: genes list with genes only in common with human
    """
    print "Filtering out genes in non-Human",
    print "species that do not exist in Humans...\n"

    # get set of common genes among all species
    otherSpecies = ["Mouse", "Chimp"]
    
    humanSet = set([]) # empty human set
    for gene in genes["Human"]: # go through human genes
        humanSet.add(gene[1].upper()) # add gene to set
    
    # get genes that in common with human pathway
    for org in otherSpecies:
        temp = set([]) # create empty set
   
        # NOTE: write commands to keep genes removed
 
        # get genes in species that are common with humans
        for gene in genes[org]: # loop through genes
            temp.add(gene[1].upper()) # add gene to a set
        genesInHuman = humanSet.intersection(temp) # find genes in human
        
        # keep genes that are common with humans
        for gene in genes[org]:
            if gene[1].upper() not in genesInHuman:
                print "Removing " + gene[1] + " from " + org
                genes[org].remove(gene) # remove gene from list

    return genes
    
def save_sequences(allAccession):
    """
    INPUT: all accession IDs for each species' genes
    OUTPUT: sequence files for each gene with sequence for each species a 
        gene exists

    The sequences will be put into the sequenceAnalysis directory that was
        created.
    """
    for gene in allAccession.keys(): # loop through genes
        print "We're going to find sequences for " + gene
        temp = [] # list to put accession numbers in
        for org in allAccession[gene].keys(): # loop through species
            temp.append(allAccession[gene][org]) # add accession to temp list
            print "Got sequence in " + org

        # only keep genes that have sequences for all three species
        if len(temp) != 3:
            print "There seems to be <3 sequences for analysis."
            print "Thus we will skip using the " + gene + " gene\n"
            continue

        handle = Entrez.efetch(db="nuccore", id=','.join(temp), 
                           rettype="fasta",retmode="text")
        print "Sequences for " + gene + " successfully obtained!"
        fasta_records = handle.read()
        handle.close()

        os.makedirs("./sequenceAnalysis/"+gene) # make folder for gene
        directory = "./sequenceAnalysis/" + gene + "/"
        fileName =  gene + ".fasta" # sequences
        directoryFile = directory + fileName # combine to get complete    

        fh = open(directoryFile, "w") # open file to put sequences
        fh.write(fasta_records) # write sequences to file
        print "Sequences written into ./sequenceAnalysis/" + gene + "\n"
        fh.close() # close file

###################
### GET PATHWAY ###
###################
"""
Take in argument from command line as target pathway
OUTPUT: pathway
This output has the name of the pathway that will be used in the analysis
"""

try:
    import sys
except ImportError:
    pass

pathway = "\"" + sys.argv[1] + "\"" # put pathway string together
print "The " + pathway + " pathway will be used for this analysis.\n"

#############################
### EXTRACT PATHWAY GENES ###
#############################
"""
Extract gene list for each species

The species that will be analyzed are:
- Homo sapiens (Humans)
- Mus musculus (Mouse)
- Pan troglodytes (Chimpanzee)

OUTPUT: genes
This output is a dictionary of lists where the key is the organism/species
(e.g. "Human") and the values are lists of lists of genes formatted like so:
[ID_NUMBER, GENE_NAME, GENE_DESCRIPTION]

The last part of this section will only keep genes from others species (Chimps
and Mice) if they are in the set of genes in Humans.
"""

try:
    import re
except ImportError:
    pass

try:
    import urllib
except ImportError:
    pass

# species legend
species = {"Human":"hsa",
           "Mouse":"mmu",
           "Chimp":"ptr"}

# dictionary for each species
genes = {}

# loop through species
for org in species.keys():
    orgList = parse_kegg_html(pathway, species[org])
    print "KEGG website was extracted for a " + org
    genes[org] = get_genes(orgList)
    print "KEGG website was parsed for " + org + " genes.\n"


# keep genes in other species that are common with humans
genes = keep_genes_common_with_humans(genes)
print "\nThe remaining number of genes from each species is:"
for org in genes.keys():
    print org + " has " + str(len(genes[org])) + " number of genes"
print "\nFinished filtering out only Human genes.\n"

#############################
### OBTAIN mRNA ACCESSION ###
#############################
"""
Get mRNA accession numbers for each gene
OUTPUT: allAccession

The output is a dictionary of dictionaries where the key is a gene and the
values are a dictionary with their key as the organism (e.g. "Human") and 
the value is the accession number for the gene from that organism 
KEY: organism
VALUE: accession
"""
print "Beginning mRNA accession requests...\n"

try:
    from Bio import Entrez
except ImportError:
    pass

try:
    import xml.etree.ElementTree as et
except ImportError:
    pass

# provide email address
email = "leunge@ohsu.edu"
Entrez.email = email # add email to object

# dictionary of dictionaries for accession IDs for each species
allAccession = {}

# loop through all species to create dictionary
# KEYS:gene, VALUE:dictionary with accession for diff species
for org in genes.keys():
    numAccess = len(genes[org]) # total num of genes to get
    i = 1 # index to keep track of

    # loop through genes
    for gene in genes[org]:
        print str(i)+"out of "+str(numAccess)+" genes in "+org
        target = gene[0] # focus on ID
        geneName = gene[1].upper() # focus on uppercase name
        print "Searching for " + geneName + " in a " + org + "..."
        if geneName not in allAccession.keys(): # if first instance
            allAccession[geneName] = {} # make dictionary
            allAccession[geneName][org] = get_accession(target)        
        else: # if this isn't first instance
            allAccession[geneName][org] = get_accession(target)
        i += 1 # increment index

print "All accession IDs for all sequences in mind have been fetched."

##########################
### CREATE DIRECTORIES ###
##########################
"""
Create directories for analysis
OUTPUT: None
The output here is to create the folder for sequence analysis
"""

try:
    import os
except ImportError:
    pass

# create folder to put sequences in if it doesn't exist
if not os.path.exists("sequenceAnalysis"):
    os.makedirs("sequenceAnalysis")

print "sequenceAnalysis directory created to save sequences and analyses\n"

#############################
### OBTAIN mRNA SEQUENCES ###
#############################
"""
Loop through accession numbers to get sequences
"""

save_sequences(allAccession)
print "Successfully saved all sequences that we want for analysis.\n"

#########################
### CLUSTAL ALIGNMENT ###
#########################

try:
    from Bio.Align.Applications import ClustalwCommandline
except ImportError:
    pass

# change path for clustalw as necessary
# clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe" # windows
# clustalw_app = r"/Applications/clustalw-2.1-macosx/clustalw2" # macintosh
# assert os.path.isfile(clustalw_app), "Clustal W executable missing"
# clustalw_cline = ClustalwCommandline(clustalw_app, infile="")
# stdout, stderr = clustalw_cline()
