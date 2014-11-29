bmi656-research-project
=======================

This repository contains the files for the final project in BMI 656.

# odds.py

This script will perform an odds ratio calculation to identify pathways
containing a larger number of differentially expressed genes than would be
expected by chance.

This script will assume you have the following text files in the current
directory:

+ H5N1_VN1203_DE_Probes.txt
+ H5N1_VN1203_UNIVERSE_Probes.txt
+ KEGG_Pathway_Genes.txt

The script will parse through each of them and find KEGG pathways which have an
odds ratio greater than 1.5. A KEGG pathway will be chose out of the ones with
an odds ratio greater than 1.5 to study further in the next section.

# crossSpecies.py
