#Resource conservation manifests in the genetic code
by Liat Shenhav and David Zeevi

DOI: http://biorxiv.org/lookup/doi/10.1101/790345

Correspondence to: dzeevi@rockefeller.edu

##Overview

**What this is.** This is the computational pipeline used to process and analyze metagenomic data used in this manuscript. By downloading the code and following the instructions below, you will be able to replicate the results reported in this work. If you are proficient in python, you will also be able to manipulate this code for a SNP-level analysis of marine microbes from metagenomic samples. 

**What this is not.** This is not software. Please do not expect to press 'GO' and get all the analyses printed out on the screen. Running this requires some know how and configuration work. It also requires significant computing power. Running the data processing and analyses below would take more than 100,000 CPU hours, a significant amount of RAM and about 35TB of space. So please do not run on your home PC unless you have 11 years to spare.


##Instructions for replication
###Get the code
`git clone https://github.com/zeevilab/resource-conservation.git` 

###Install required software, packages and dependencies
* Python packages: 
	* pandas
	* scipy
	* dill
	* pysam
	* biopython 
* External programs: 
	* bowtie2
	* samtools
	* bcftools
	* eggNOG-mapper
* Build cython files in directory Data/Mapping/cy_ext_seq by running setup.py (this may raise an error if your system is missing a c++ compiler - do not ignore it). 

###Get the data


###Get the metadata


###Get the OM-RGC database, index it and run eggNOG-mapper on it


###Edit the config file



###Configure your HPC system



###Run data creation pipelines
####MappingPipe.py
####CallingPipe.py
####SNPPipe.py



###Run the analysis pipelines
