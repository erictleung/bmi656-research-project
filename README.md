bmi656-research-project
=======================

This repository contains the files for the final project in BMI 656.

**Table of Contents**

- [Prerequisites](#prerequisites)
- [`odds.py`](#odds.py)
- [`crossSpecies.py`](#crossspecies.py)
- [Helpful Resources](#helpful-resources)

Prerequisites
--------------

- Python 3
- Python 3 Modules
    - `os`
    - `re`
    - `urllib`
    - `Bio` (BioPython)
    - `xml`
    - `matplotlib`
- Clustal W (http://www.clustal.org/clustal2/)

odds.py
-------

This script will perform an odds ratio calculation to identify pathways
containing a larger number of differentially expressed genes than would be
expected by chance.

This script will assume you have the following text files in the current
directory:

- `H5N1_VN1203_DE_Probes.txt`
- `H5N1_VN1203_UNIVERSE_Probes.txt`
- `KEGG_Pathway_Genes.txt`

The script will parse through each of them and find KEGG pathways which have an
odds ratio greater than 1.5. A KEGG pathway will be chose out of the ones with
an odds ratio greater than 1.5 to study further in the next section.

crossSpecies.py
---------------

Deliverables
------------

- (**60 points**) Write Part I and Part II as separate python programs that are
  both called from a linux bash script. This script should check for the
  existence of relevant files and check that Part I has successfully completed
  before running Part II. Turn in both python programs and the bash script.
- (**10 points**) Screen shot of pathway with highlighted differentially
  expressed genes
- (**10 points**) Boxplot labeled with p-value of **Mann-Whitney statistic**
- (**20 points**) 1-2 page write-up summarizing your findings. This should be a
  word document that includes at least 2 figures corresponding to your pathway
  and boxplot. Discuss the limitations of the study and any obstacles/problems
  you experienced. Discuss the relevance of the pathway to H5N1 infection.
  Discuss why the conservation of the affected pathway might be important for
  the study of H5N1 (Hint: think about our use of model organisms to study
  infectious diseases).

Helpful Resources
-----------------

- [Biopython Tutorial and Cookbook][biopython]

[biopython]: http://biopython.org/DIST/docs/tutorial/Tutorial.html
