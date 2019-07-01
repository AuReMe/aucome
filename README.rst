.. image:: https://img.shields.io/pypi/v/aucome.svg
	:target: https://pypi.python.org/pypi/aucome

AuCoMe: Automatic Comparison of Metabolism
==========================================

Workflow to reconstruct multiple metabolic networks in order to compare them.

.. contents:: Table of contents
   :backlinks: top
   :local:


Installation
------------

Dependencies
~~~~~~~~~~~~

This tool needs:

	- `Orthofinder <https://github.com/davidemms/OrthoFinder>`__ (which needs `Diamond <https://github.com/bbuchfink/diamond>`__ and `FastME <https://gite.lirmm.fr/atgc/FastME/>`__)

	- `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ (which needs Blast)

And some python package:

	- `biopython <https://github.com/biopython/biopython>`__

	- `docopt <https://github.com/docopt/docopt>`__

	- `eventlet <https://github.com/eventlet/eventlet>`__

	- `intervene <https://github.com/asntech/intervene>`__ (which needs R package `UpSetR <https://github.com/hms-dbmi/UpSetR>`__)

	- `lxml <https://github.com/lxml/lxml>`__

	- `matplotlib <https://github.com/matplotlib/matplotlib>`__

	- `mpwt <https://github.com/AuReMe/mpwt>`__

	- `pandas <https://github.com/pandas-dev/pandas>`__

	- `padmet <https://github.com/AuReMe/padmet>`__

	- `padmet-utils <https://github.com/AuReMe/padmet-utils>`__

	- `requests <https://github.com/kennethreitz/requests>`__

	- `seaborn <https://github.com/mwaskom/seaborn>`__


To run annotation based reconstruction, you need to install Pathway Tools. This tool is available at the `Pathway Tools website <http://bioinformatics.ai.sri.com/ptools/>`__. A command in the package install the tools:

.. code:: sh

        aucome --installPWT=path/to/pathway/tools/installer

You can also provide an option to this commande: --ptools=ptools_path

This option let you choose the path where the ptools-local folder will be installed. PGDBs created by Pathway-Tools are stored in this folder.


Docker
~~~~~~

From git repository:

.. code:: sh

	git clone https://github.com/AuReMe/aucome.git

	cd aucome

	docker build -t "my_image".


pip
~~~

If you have installed all the dependencies, you can just install acuome with:

.. code:: sh

	pip install aucome

Usage
-----

Initialization
~~~~~~~~~~~~~~

You have to create the working folder for AuCoMe, with the --init argument:

.. code:: sh

    aucome --init=run_ID [-v]

This command will create a folder name "run_ID" inside the working folder. In this "run_ID" folder, the command will create all the folders used during the analysis.

.. code-block:: text

	run_ID
	├── analysis
		├── group_template.tsv
		├──
	├── annotation_based
		├── PADMETs
			├──
		├── PGDBs
			├──
		├── SBMLs
			├──
	├── config.txt
	├── logs
		├──
	├── model_organisms
		├──
	├── networks
		├── PADMETs
			├──
		├── SBMLs
			├──
	├── orthology_based
		├── Orthofinder_WD
			├──
	├── studied_organisms
		├──

**analysis** will store the result of padmet analysis.

**annotation_based** contains three sub-folders. The folder PGDBs will contain all the results from Pathway-Tools (in dat format). These results will be stored in padmet and sbml inside PADMETs and SBMLs.

**config.txt** contains numerous paths used by the script.

**model_organisms** contains the model organisms you want to use for the orthology. In this folder you put a new folder with the name of the species and in this folder you put the proteome and the sbml of the metabolic network of your species. Proteome and metabolic network names must be the same than the name of the folder.

.. code-block:: text

	├── model_organisms
		├── A_thaliana
			├── A_thaliana.fasta
			├── A_thaliana.sbml

**networks** will contain all the metabolic network created by aucome in padmet format.

**orthology_based** contains one folder Orthofinder_WD. This folder will contain all the run of Orthofinder.

**studied_organisms**: you put all the species that you want to studies in this folder. For each species you create a folder and in this folder you put the genbank file of this species. Like for model_organisms, file and folder must have the same name. And the genbank file must end with a '.gbk'.

.. code-block:: text

	├── studied_organisms
		├── species_1
			├── species_1.gbk
		├── species_2
			├── species_2.gbk


Once you have put your species in the studied_organisms folder and the model in model_organisms, a check must be done on the data using:

Check command
~~~~~~~~~~~~~

.. code:: sh

    aucome check --run=run_ID [--cpu=INT] [-v]

This command will check if there is no character that will cause trouble. It will also create the proteome fasta file from the genbank.

Also, this command will fill the 'all' row of analysis/group_template.tsv, with all the species from the studied_organisms folder.

And for the annotation_based folder, if PGDBs contains folder, it will create the padmet and the sbml corresponding to these draft in PADMETs and SBMLs.

Renconstruction command
~~~~~~~~~~~~~~~~~~~~~~~

A run of Pathway-Tools can be launched using the command:

.. code:: sh

    aucome reconstruction --run=run_ID [--cpu=INT] [-v]

.. code-block:: text

	├── annotation_based
		├── PADMETs
			├── output_pathwaytools_species_1.padmet
			├── output_pathwaytools_species_2.padmet
		├── PGDBs
			├── species_1
				├── PGDB dat files
				├── ...
			├── species_2
				├── PGDB dat files
				├── ...
		├── SBMLs
			├── output_pathwaytools_species_1.sbml
			├── output_pathwaytools_species_2.sbml
	├── logs
		├── log_error.txt
		├── resume_inference.tsv

Using the package mpwt, it will create the input file for Pathway-Tools inside studied_organisms and if there is no error, it will create for each species inside this folder a folder inside PGDBs containing all the dat files ofthe draft metabolic network.

Orthology command
~~~~~~~~~~~~~~~~~

Orthofinder can be launched using:

.. code:: sh

	aucome orthology --run=run_ID [-S=STR] [--orthogroups] [--cpu=INT] [-v]

.. code-block:: text

	├── orthology_based
		├── Orthofinder_WD
			├── species_1
				├── output_orthofinder_from_species_2.sbml
			├── species_2
				├── output_orthofinder_from_species_1.sbml
			├── Orthofinder_WD
				├── species_1.faa
				├── species_2.faa
				├── OrthoFinder
					├── Results_MonthDay
						├── Orthogroups
						├── Orthologues
						├── ..

Then the proteome from the studied organisms and from the models will be moved to the Orthofinder_WD folder and orthofinder will be launch on them. Orthofinder result will be in this folder and in orthology_based, there will be all the metabolic network reconstructed from orthology.

Draft command
~~~~~~~~~~~~~

Then you can merge all the metabolic network with:

.. code:: sh

    aucome draft --run=run_ID [--cpu=INT] [-v]

.. code-block:: text

	├── networks
		├── PADMETs
			├── species_1.padmet
			├── species_2.padmet
		├── SBMLs
			├── species_1.sbml
			├── species_2.sbml

This will output the result inside the networks folder.

Workflow command
~~~~~~~~~~~~~~~~

You can launch the all workflow with the command:

.. code:: sh

    aucome workflow --run=ID [-S=STR] [--orthogroups] [--cpu=INT] [-v]

Analysis command
~~~~~~~~~~~~~~~~

You can launch group analysis with the command:

.. code:: sh

    aucome analysis --run=ID [--cpu=INT] [-v]

You must write the groups of species that you want to analyze in the analysis/group_template.tsv file:
The first line of the file contains 'all' (it will launch the analysis on all the species).

When you create the repository with --init, the file will only contain 'all' row:

+--------------+------------+-------------+--------------+--------------+
|   all        |            |             |              |              |
+--------------+------------+-------------+--------------+--------------+

After the check (with check or workflow command), it will add all the species that you have in your studied_organisms folder:

+--------------+------------+-------------+--------------+--------------+
|   all        | species_1  | species_2   | species_3    | species_4    |
+--------------+------------+-------------+--------------+--------------+

Then you can create a new row to add another group. The name of the group is in the first column. Then for each species you add a column with the species name.
You must at least give 2 species.

Example:

+--------------+------------+-------------+--------------+--------------+
|   all        |species_1   | species_2   | species_3    | species_4    |
+--------------+------------+-------------+--------------+--------------+
|   group_1    | species_1  | species_2   |              |              |
+--------------+------------+-------------+--------------+--------------+
|   group_2    | species_1  | species_2   | species_4    |              |
+--------------+------------+-------------+--------------+--------------+

This script will create one folder for each group:

.. code-block:: text

	├── analysis
		├── group_template.tsv
		├── all
			├──
		├── group_1
			├──
		├── group_2
			├──

Compare command
~~~~~~~~~~~~~~~~

You can launch group analysis with the command:

.. code:: sh

    aucome analysis --run=ID [--cpu=INT] [-v]

This script will read the group_template.tsv file and create a folder containing an upset graph comparing the group that you selected:

.. code-block:: text

	├── analysis
		├── group_template.tsv
		├── upgset_graph
			├── genes.csv
			├── Intervene_upset.R
			├── Intervene_upset.svg
			├── Intervene_upset_combinations.txt
			├── metabolites.csv
			├── pathways.csv
			├── reactions.csv
			├── tmp_data
				├──