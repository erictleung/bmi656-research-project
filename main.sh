#!/bin/bash

# Name:       Eric Leung
# Assignment: Research Project
# Date:       December 8th, 2014
# Script:     main.sh
# Desription: This is the main command to 

# perform OR on pathways and choose pathway
pathway=$(python odds.py)

# examine cross-species conservation of pathway genes
python crossSpecies.py $pathway
