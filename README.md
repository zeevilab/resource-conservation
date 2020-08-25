# Resource conservation manifests in the genetic code
by Liat Shenhav and David Zeevi

DOI: http://biorxiv.org/lookup/doi/10.1101/790345

Correspondence to: lshenhav@rockefeller.edu, dzeevi@rockefeller.edu

## Overview

**What this is.** This is the computational pipeline used to process and analyze metagenomic data used in this manuscript. By downloading the code and following the instructions below, you will be able to replicate the results reported in this work. If you are proficient in python, you will also be able to manipulate this code for a SNP-level analysis of marine microbes from metagenomic samples. 

## Instructions for replication
### Get the code
`git clone https://github.com/zeevilab/resource-conservation.git` 

### Install required software, packages and dependencies
* This pipeline was run with Python 3.6.9 (Anaconda). eggNOG mapper was run using Python 2.7.
* Python packages: 
	* numpy (1.17.2)
	* pandas (0.25.1)
	* scipy (1.3.1)
	* dill (0.3.1.1)
	* pysam (0.15.3)
	* biopython (1.74)
	* xarray (0.10.7)
* External programs: 
	* bowtie2 (2.3.2)
	* samtools (1.8)
	* bcftools (1.6)
	* eggNOG-mapper
* Build cython files in directory Data/Mapping/cy_ext_seq by running setup.py (this may raise an error if your system is missing a C++ compiler - do not ignore it). 

### Get the data
Download metagenomic data (fastq files) for Tara oceans (prokaryote size specification), bioGEOTRACES and HOT/BATS from ENA. Also save the sample sheet for these samples. Accession numbers:
* ENA:PRJEB1787 (TARA oceans prokaryotic fraction) 
* ENA:PRJNA385854 (bioGEOTRACES) 
* ENA:PRJNA385855 (HOT/BATS)

### Get the metadata
Download metadata for Tara oceans, GEOTRACES, HOT and BATS into separate folders and link to them in the code (instructions to follow):
* HOT data from station ALOHA (http://aloha-data.soest.hawaii.edu/repository/entry/show/Home/HOT+Data/Station+02+(ALOHA)): 
	* Download CTD and water data into two separate folders.
* BATS data from http://bats.bios.edu/bats-data/:
	* Download bottle data (bats_bottle.txt).
	* Download all ASCII records into a separate folder.
* GEOTRACES data from https://www.bodc.ac.uk/geotraces/data/idp2017/:
	* Download GEOTRACES IDP2017 v2 discrete sample data in ASCII format.
	* Download GEOTRACES IDP2017 v1 CTD sensor data in ASCII format.
* TARA oceans data from PANGAEA:
	* Samples context sequencing: https://doi.pangaea.de/10.1594/PANGAEA.875581
	* Carbonate chemistry: https://doi.pangaea.de/10.1594/PANGAEA.875567
	* Nutrients: https://doi.pangaea.de/10.1594/PANGAEA.875575
	* HPLC: https://doi.pangaea.de/10.1594/PANGAEA.875569
	* Station registry: https://doi.pangaea.de/10.1594/PANGAEA.842227

### Get the OM-RGC database and index it
* Get the OM-RGC database here: http://ocean-microbiome.embl.de/companion.html as a tab-delimited file.
* Create a FASTA file from the tab delimited file containing the samples, where the id of each sequence is the OM-RGC_ID (column 1) and the sequence is the sequence (column 13), as specified here: http://ocean-microbiome.embl.de/data/OM-RGC_Readme.release.txt.
* Index the FASTA file using bowtie2_build

### Edit the config file
The configuration file config.py holds all dataset locations, metadata locations, locations for temporary files and parameters used. Edit config.py as directed in comments within the file to point to your downloaded resources above. *It is very important to take care linking to the correct resources as all the following steps depend on this*.

### Configure your HPC system
Running the data creation pipelines can be a lengthy process. It is advisable to configurate it to run on a multi-CPU machine with sufficient memory. We ran it on a computational cluster with 56 nodes, each with 24 CPUs and 256GB RAM, and an additional 'fat' node with 64 CPUs and 3TB RAM.
As different HPCs run different job scheduling software, we have noted in our data creation pipelines the labor intensive positions in the code where a *for* loop should be replaced with adequate HPC job scheduling to significantly cut running time. Follow comments in Pipe files listed below.

### Run data creation pipelines
#### MetadataPipe.py
The metadata and annotation pipeline prepares all the necessary data for running the following pipelines. It has two main parts parts:
1. Concentrating all the relevant metadata in one place:
	1. It reads the sample sheet saved from ENA and creates a reference for all files to be used in the pipeline.
	1. It takes the relevant environmental measurements from the data downloaded from Tara, GEOTRACES, ALOHA and BATS.
1. Annotating the OM-RGC database with eggnog-mapper.

#### MappingPipe.py
The mapping pipeline aligns reads from all sources to the OM-RGC database. It then runs ICRA on it (Zeevi et al., Nature 2019) to correct read assignment and updates BAM files resulting from the mapping with the probability of alignment for each read to each OM-RGC gene.

#### CallingPipe.py
The calling pipeline sorts and filters each BAM file, and then splits it into chunks according to genes to which reads map. It then runs bcftools on each of the chunks to call variants across all samples.

#### SNPPipe.py
The SNP pipeline filters variants according to filtering thresholds defined in the config file. For each gene, it calculates population genomic markers (pN/pS, pi) and then unites genes according to KEGG and eggNOG orthologous groups. The output of this stage are matrices for each of the markers and each of the OG databases (KEGG/eggNOG).

### Run the analysis pipelines
#### Linear mixed models (LMM)
Lorem ipsum
#### Genetic code analyses
The python file CodeAnalysisPipe.py in directory Analysis/GeneticCode groups all the analyses performed
regarding the genetic code. These should be run only after running all of the data creation pipelines.
In addition, a configuration of the HPC system is required to run these analyses.