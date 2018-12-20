---
title:  compare_metabo - Documentation
author: Meziane AITE & Arnaud BELCOUR
date: 2018-07-13
version: 0.4
---

################################################################################

## Description

Workflow to reconstruct multiple metabolic networks in order to compare them.

## Installation

From git repository:

	git clone https://gitlab.inria.fr/DYLISS/compare_metabo.git

	cd compare_metabo

	docker build .

To run annotation based reconstruction, you need to install Pathway-Tools. This tool is available at the [Pathway-Tools website](http://bioinformatics.ai.sri.com/ptools/). A command in the package install the tools:

        compare --installPWT=path/to/pathway/tools/installer

You can also provide an option to this commande: --ptools=ptools_path

This option let you choose the path where the ptools-local folder will be installed. PGDBs created by Pathway-Tools are stored in this folder.


## Architecture

        .
        ├── README.md
        ├── compare.py
        ├── Dockerfile
        ├── update_version.py
        ├── release.txt

## Usage

First you must choose the working folder for the tool. By default, it is "/shared". This can be achieve by the commande:

        compare --setWorkingFolder=DIR

Once the working folder set, you can create a run folder of compare using the command init:

        compare --init=run_ID [-v]

This command will create a folder name "run_ID" inside the working folder. In this "run_ID" folder, the command will create all the folders used during the analysis.

        run_ID
        ├── analysis
             ├──
        ├── annotation_based
             ├── PADMETs
                 ├──
             ├── PGDBs
                 ├──
             ├── SBMLs
                 ├──
        ├── config.txt
        ├── model_organisms
             ├──
        ├── networks
             ├──
        ├── orthology_based
             ├── Orthofinder_WD
                 ├──
        ├── studied_organisms
             ├──

analysis will store the result of padmet analysis.

annotation_based contains three sub-folders. The folder PGDBs will contain all the results from Pathway-Tools (in dat format). These results will be stored in padmet and sbml inside PADMETs and SBMLs.

config.txt contains numerous paths used by the script.

model_organisms contains the model organisms you want to use for the orthology. In this folder you put a new folder with the name of the species and in this folder you put the proteome and the sbml of the metabolic network of your species. The name of the proteome and of the metabolic network must the same than the name of the folder.

        ├── model_organisms
             ├── A_thaliana
                 ├── A_thaliana.fasta
                 ├── A_thaliana.sbml

networks will contain all the metabolic network created by compare in padmet format.

orthology_based contains one folder Orthofinder_WD. This folder will contain all the run of Orthofinder.

studied_organisms: you put all the species that you want to studies in this folder. For each species you create a folder and in this folder you put the genbank file of this species. Like for model_organisms, file and folder must have the same name. And the genbank file must end with a '.gbk'.

        ├── studied_organisms
             ├── species_1
                 ├── species_1.gbk
             ├── species_2
                 ├── species_2.gbk


Once you have put your species in the studied_organisms folder and teh model in model_organisms, a check must be done on the data using:

        compare --run=run_ID -c [-v]

This command will check if there is no character that will make some scritp crashed later in the analysis. It will also create the proteome fasta file from the genbank.

And for the annotation_based folder, if PGDBs contains folder, it will create the padmet and the sbml corresponding to these draft in PADMETs and SBMLs.

A run of Pathway-Tools can be launched using the command:

        compare --run=run_ID -p [-v]

Using the package mpwt, it will create the input file for Pathway-Tools inside studied_organisms and if there is no error, it will create for each species inside this folder a folder inside PGDBs containing all the dat files ofthe draft metabolic network.

Orthofinder can be launched using:

        compare --run=run_ID -o [-v]

Then the proteome from the studied organisms and from the models will be moved to the Orthofinder_WD folder and orthofinder will be launch on them. Orthofinder result will be in this folder and in orthology_based, there will be all the metabolic network reconstructed from orthology.

Then you can merge all the metabolicnetwork with:

        compare --run=run_ID -d [-v]

This will output the result inside the networks folder.
