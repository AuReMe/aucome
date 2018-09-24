---
title:  compare_metabo - Documentation
author: Meziane AITE & Arnaud BELCOUR
date: 2018-07-13
version: 0.1
---

################################################################################

## Description

Workflow to reconstruct multiple metabolic networks in order to compare them.

## Installation

From git repository:

	git clone https://gitlab.inria.fr/DYLISS/compare_metabo.git

	cd compare_metabo

	docker build .

To run annotation based reconstruction, you need to install Pathway-Tools. This tool is available at the [Pathway-Tools website](http://bioinformatics.ai.sri.com/ptools/). Then run the command:

        compare --installPWT=path/to/pathway/tools/installer

## Architecture

        .
        ├── README.md
        ├── compare.py
        ├── ptools

