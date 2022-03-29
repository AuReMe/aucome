.. image:: https://img.shields.io/pypi/v/aucome.svg
	:target: https://pypi.python.org/pypi/aucome

AuCoMe: Automatic Comparison of Metabolism
==========================================

**WORK IN PROGRESS** Workflow to reconstruct multiple metabolic networks in order to compare them.

.. contents:: Table of contents
   :backlinks: top
   :local:


Installation
------------

Dependencies
~~~~~~~~~~~~

These tools are needed:

	- `Exonerate <https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate>`__

	- `Orthofinder <https://github.com/davidemms/OrthoFinder>`__ (which needs `Diamond <https://github.com/bbuchfink/diamond>`__, `FastME <https://gite.lirmm.fr/atgc/FastME/>`__, and `MMseqs2 <https://github.com/soedinglab/MMseqs2/>`__)

	- `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ (which needs `Blast <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__)

	- `R <https://cran.r-project.org/>`__

And some python packages:

	- `biopython <https://github.com/biopython/biopython>`__
	
	- `cobra <https://github.com/opencobra/cobrapy>`__ 
	
	- `cython <https://github.com/cython/cython>`__

	- `docopt <https://github.com/docopt/docopt>`__

	- `eventlet <https://github.com/eventlet/eventlet>`__

	- `lxml <https://github.com/lxml/lxml>`__

	- `matplotlib <https://github.com/matplotlib/matplotlib>`__

	- `mpwt <https://github.com/AuReMe/mpwt>`__

	- `numpy <https://github.com/numpy/numpy>`__

	- `padmet <https://github.com/AuReMe/padmet>`__

	- `pandas <https://github.com/pandas-dev/pandas>`__

	- `pybind11 <https://github.com/pybind/pybind11>`__ 

	- `pyparsing <https://github.com/pyparsing/pyparsing>`__

	- `pythran <https://github.com/serge-sans-paille/pythran>`__ 

	- `requests <https://github.com/kennethreitz/requests>`__

	- `rpy2 <https://github.com/rpy2/rpy2>`__

	- `scipy <https://github.com/scipy/scipy>`__

	- `seaborn <https://github.com/mwaskom/seaborn>`__

	- `setuptools <https://github.com/pypa/setuptools>`__

	- `supervenn <https://github.com/gecko984/supervenn>`__

	- `tzlocal <https://github.com/regebro/tzlocal>`__

Installation of Pathway Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run annotation based reconstruction, you need to install Pathway Tools. This tool is 
available at the `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ website. A 
command in the package install the tools:

.. code:: sh

        aucome --installPWT=path/to/pathway/tools/installer
	source ~/.bashrc

You can also provide an option to this commande: --ptools=ptools_path


This option let you choose the path where the ptools-local folder will be installed. PGDBs 
created by `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ are stored in this 
folder.


Getting the MetaCyc PADMET file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You also should install the MetaCyc_XX.X.padmet (the version number of 
`MetaCyc <https://metacyc.org/>`__  is replaced with XX.X), and then you should update your 
config.txt files for each study. This is the way to 
getting a MetaCyc_XX.padmet file: Firstly, download the flat files of 
`MetaCyc <https://metacyc.org/>`__ in DAT format at the
`https://biocyc.org/download.shtml <https://biocyc.org/download.shtml>`__ webpage. Secondly, 
put all the downloaded DAT files in a directory (it is named FLAT_DIR here). Thirdly run this 
command:

.. code:: sh

	padmet pgdb_to_padmet --pgdb=FLAT_DIR --output=metacyc_XX.X.padmet --version=XX.X --db=Metacyc -v


Docker
~~~~~~

From git repository:

.. code:: sh

	git clone https://github.com/AuReMe/aucome.git

	cd aucome

	docker build -t "my_image".


Singularity
~~~~~~~~~~~

You need to have a pathway tools installer on the same path as the recipe.

From git repository:

.. code:: sh

	sudo singularity build aucome.sif Singularity

If you have the issue:

.. code:: sh

	FATAL:   While performing build: while creating squashfs: create command failed: exit status 1: Write failed because No space left on device
	FATAL ERROR: Failed to write to output filesystem

It is because Singularity has not enough space in its temporary folder due to the size of the
tools needed by aucome. You can modify manually this path using the ``SINGULARITY_TMPDIR`` 
variable (the temporary folder must exist), for example:

.. code:: sh

	sudo SINGULARITY_TMPDIR=/home/user/tmp_folder singularity build  aucome.sif Singularity

Then you can run the container with command like:

.. code:: sh

	singularity run  aucome.sif aucome workflow --run data  --filtering --cpu 10

But using only these commands can produce errors due to the compartmentalization of singularity.
So it is better to use the ``-c`` to avoid sharing filesystem with host.
And the ``-B`` allows to give a shared folder between the host and the singularity container 
so Singularity can also access to the data in the host.

.. code:: sh

	singularity run -c -H /path/outside/singularity/to/shared:/path/in/singularity/container aucome.sif aucome workflow --run /path/in/singularity/container/data  --filtering --cpu 10


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

This command will create a folder name "run_ID" inside the working folder. In this "run_ID"
folder, the command will create all the folders used during the analysis.

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
	├── networks
		├── PADMETs
			├──
		├── SBMLs
			├──
	├── orthology_based
		├── 0_Orthofinder_WD
			├── OrthoFinder
		├── 1_sbml_orthology
		├── 2_padmet_orthology
		├── 3_padmet_filtered
	├── structural_check
		├── 0_specifics_reactions
		├── 1_blast_results
			├── analysis
			├── tmp
		├── 2_reactions_to_add
		├── 3_PADMETs
	├── studied_organisms
		├──

**analysis** will store the various analysis of the 
`PADMET <https://padmet.readthedocs.io/en/latest/tutorial.html#padmet-format>`__ files which 
are in the networks folder.

**annotation_based** includes three subfolders. The PGDBs folder will contain all the results 
from `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ (in DAT format). These results
will also be stored in `PADMET <https://padmet.readthedocs.io/en/latest/tutorial.html#padmet-format>`__ 
and `SBML <https://sbml.org/documents/specifications/>`__ files inside PADMETs and SBMLs.

**config.txt** contains numerous paths used by the script: paths to programs, directories and 
databases. It also inclues the `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ 
and `MetaCyc <https://metacyc.org/>`__  versions. 

**networks** will contain one metabolic network per studied organism, created thanks to AuCoMe,
in `PADMET <https://padmet.readthedocs.io/en/latest/tutorial.html#padmet-format>`__ and 
`SBML <https://sbml.org/documents/specifications/>`__ formats that are stored into two
directories (PADMETs and SBMLs). It also includes the panmetabolism of all the studied 
organisms in `PADMET <https://padmet.readthedocs.io/en/latest/tutorial.html#padmet-format>`__
and `SBML <https://sbml.org/documents/specifications/>`__ format. 

**orthology_based** contains four subfolders. Firstly the 0_Orthofinder_WD directory folder 
will include all the run of `Orthofinder <https://github.com/davidemms/OrthoFinder>`__. 
Secondly, the 1_sbml_orthology folder will contain one subdirectory per studied organims, and 
each subfolders include `SBML <https://sbml.org/documents/specifications/>`__  files with the
orthogroups of other species that `OrthoFinder <https://github.com/davidemms/OrthoFinder>`__ 
found. Thirdly, the 2_padmet_orthology directory will contain the 
`PADMET <https://padmet.readthedocs.io/en/latest/tutorial.html#padmet-format>`__ files created 
with the orthology step. Fourthly, the 3_padmet_filtered folder will contain 
`PADMET <https://padmet.readthedocs.io/en/latest/tutorial.html#padmet-format>`__ files created
thanks to the orthology step, but in this subfolder only the robust reactions are kept in 
these `PADMET <https://padmet.readthedocs.io/en/latest/tutorial.html#padmet-format>`__ files.  

**structral_check** relies on the search on the genomes for missing Gene-Proteins-Reactions 
associations. All the metabolic networks previously created are be pairwise compared. If one 
metabolic network has a Gene-Protein-Reaction association that another one has not, a genomic 
search will be performed between both genomes corresponding with the both metabolic networks.
Gene-Protein-Reaction associated with the first metabolic network will be used to search for 
match with the genome sequence corresponding with of the second metabolic network.
It contains four subdirectories. Firstly 0_specifics_reactions folder will include numerous 
TSV files with lists of Gene-Protein-Reaction associations that are present in a metabolic 
network and that are absent in another metabolic network. Secondly, the 1_blast_results 
directory will contain the search results between genomes of studied organisms and selected 
genes in the previous TSV files. Here orther TSV files will also be created with another format. These TSV 
files will include the results of genomic search programs. 
`BlastP <https://blast.ncbi.nlm.nih.gov/>`__, `TblastN <https://blast.ncbi.nlm.nih.gov/>`__, 
and `Exonerate <https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate>`__ are 
used as genomic search programs. Thirdly the 2_reactions_to_add folder will contain a PADMET 
form with the reactions to add for each studied organisms. Fourthly, the 3_PADMETs will include
the `PADMET <https://padmet.readthedocs.io/en/latest/tutorial.html#padmet-format>`__ files 
created with the structural step.

**studied_organisms**: you put all the species that you want to study in this folder. For each 
species, you create a folder and in this folder you put the 
`GenBank <https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>`__ file of this species. Each
files and folders must have the same name. Then, the 
`GenBank <https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>`__ file must end with a 
'.gbk'.

.. code-block:: text

	├── studied_organisms
		├── species_1
			├── species_1.gbk
		├── species_2
			├── species_2.gbk


.. warning:: Remember to check the versions of `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ and `MetaCyc <https://metacyc.org/>`__ before running the check command. 

Once you have put your species in the studied_organisms folder, a check must be done on the data using:

Check command
~~~~~~~~~~~~~

.. code:: sh

    aucome check --run=ID [--cpu=INT] [-v] [--vv]

This command will check if there is no character that will cause trouble. It will also create
the proteome `FASTA <http://bioinformatics.org/annhyb/examples/seq_fasta.html>`__ file from 
the `GenBank <https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>`__. Also, this command
will fill the 'all' row of analysis/group_template.tsv, with all the species from the 
studied_organisms folder. And for the annotation_based folder, if PGDBs contains folder, it 
will create the `PADMET <https://padmet.readthedocs.io/en/latest/tutorial.html#padmet-format>`__
and the `SBML <https://sbml.org/documents/specifications/>`__ corresponding to these draft in 
PADMETs and SBMLs folders.

Renconstruction command
~~~~~~~~~~~~~~~~~~~~~~~

A run of Pathway-Tools can be launched using the command:

.. code:: sh

    aucome reconstruction --run=ID [--cpu=INT] [-v] [--vv]

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

	aucome orthology --run=ID [-S=STR] [--orthogroups] [--cpu=INT] [-v] [--vv] [--filtering] [--threshold=FLOAT]

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

Structural command
~~~~~~~~~~~~~~~~~~

To assure that no reactions are missing due to missing gene structures a genomic search is performed for all reactions appearing in one organism but not in another.

.. code:: sh

    aucome structural --run=ID [--keep-tmp] [--cpu=INT] [-v]

.. code-block:: text

	├── structural_check
		├── 0_specifics_reactions
			├── species_1_VS_species_2.tsv
			├── species_2_VS_species_1.tsv
		├── 1_blast_results
			├── analysis
				├── species_1_VS_species_2.tsv
				├── species_2_VS_species_1.tsv
			├── tmp
		├── 2_reactions_to_add
			├── species_1.tsv
			├── species_2.tsv
		├── 3_PADMETs
			├── species_1.padmet
			├── species_2.padmet


Spontaneous command
~~~~~~~~~~~~~

Then you can spontaneous all the metabolic network with:

.. code:: sh

    aucome spontaneous --run=ID [--cpu=INT] [-v] [--vv]

.. code-block:: text

	├── networks
		├── PADMETs
			├── species_1.padmet
			├── species_2.padmet
		├── panmetabolism.padmet
		├── panmetabolism.sbml
		├── SBMLs
			├── species_1.sbml
			├── species_2.sbml

This will output the result inside the networks folder.

Workflow command
~~~~~~~~~~~~~~~~

You can launch the all workflow with the command:

.. code:: sh

    aucome workflow --run=ID [-S=STR] [--orthogroups] [--keep-tmp] [--cpu=INT] [-v] [--vv] [--filtering] [--threshold=FLOAT]

Analysis command
~~~~~~~~~~~~~~~~

You can launch group analysis with the command:

.. code:: sh

    aucome analysis --run=ID [--cpu=INT] [--pvclust] [-v]

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

    aucome compare --run=ID [--cpu=INT] [-v]

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
