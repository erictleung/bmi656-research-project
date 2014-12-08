#!/bin/bash

# Name:       Eric Leung
# Assignment: Research Project
# Date:       December 8th, 2014
# Bash Ver:   3.2.53
# Script:     main.sh
# Desription: This is the main command to 
# check for files needed for scripts ahead
# Example:    sh main.sh

# perform OR on pathways and choose pathway
if [ -f "odds.py" ]; then
    echo "Script for first part of project exists."
    python odds.py
else
    echo "Script for first part doesn't exist."
fi

# examine cross-species conservation of pathway genes
if [ -f "pathway_info.csv" ] && [ -f "crossSpecies.py" ]; then
    echo "Script and CSV file for second part of project exists!\n"
    python crossSpecies.py pathway_info.csv
else
    echo "Script for second part of project doesn't exist."
    echo "Or first part of script didn't finish"
    echo "Please try again."
fi
