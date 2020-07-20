# Resource conservation manifests in the genetic code
by Liat Shenhav and David Zeevi

DOI: http://biorxiv.org/lookup/doi/10.1101/790345

Correspondence to: dzeevi@rockefeller.edu

## Overview

**What this is.** This is the computational pipeline used to process and analyze metagenomic data used in this manuscript. By downloading the code and following the instructions below, you will be able to replicate the results reported in this work. If you are proficient in python, you will also be able to manipulate this code for a SNP-level analysis of marine microbes from metagenomic samples. 

**What this is not.** This is not software. Please do not expect to press 'GO' and get all the analyses printed out on the screen. Running this requires some know how and configuration work. It also requires significant computing power. Running the data processing and analyses below would take more than 100,000 CPU hours, a significant amount of RAM and about 35TB of space. So please do not run on your home PC unless you have 11 years to spare.


## Instructions for replication
### Get the code
`git clone https://github.com/zeevilab/resource-conservation.git` 

### Install required software, packages and dependencies
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
	* Download GEOTRACES IDP2017 v2 CTD data in ASCII format.
* TARA oceans data from PANGAEA:
	* Samples context sequencing: https://doi.pangaea.de/10.1594/PANGAEA.875581
	* Carbonate chemistry: https://doi.pangaea.de/10.1594/PANGAEA.875567
	* Nutrients: https://doi.pangaea.de/10.1594/PANGAEA.875575
	* HPLC: https://doi.pangaea.de/10.1594/PANGAEA.875569
	* Station registry: https://doi.pangaea.de/10.1594/PANGAEA.842227

### Get the OM-RGC database, index it and run eggNOG-mapper on it
* Get the OM-RGC database here: http://ocean-microbiome.embl.de/companion.html as a tab-delimited file.
* Create a FASTA file from the tab delimited file containing the samples, where the id of each sequence is the OM-RGC_ID (column 1) and the sequence is the sequence (column 13), as specified here: http://ocean-microbiome.embl.de/data/OM-RGC_Readme.release.txt


### Edit the config file



### Configure your HPC system



### Run data creation pipelines
#### MappingPipe.py
#### CallingPipe.py
#### SNPPipe.py



### Run the analysis pipelines
