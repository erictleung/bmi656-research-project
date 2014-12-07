#!/bin/bash
# Name:       Eric Leung
# Assignment: Research Project
# Date:       December 8th, 2014
# Bash Ver:   3.2.53
# Script:     main.sh
# Desription: This is the main command to 
# check for files needed for scripts ahead
if [ -f "Dodds.py" ]; then
    echo "Script for first part of project exists."
else
    echo "Script for first part doesn't exist."
fi
# perform OR on pathways and choose pathway
# examine cross-species conservation of pathway genes
# python crossSpecies.py
