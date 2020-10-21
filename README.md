# Snakemake workflow: Phage_clouds_snake

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/Phage_clouds_snake.svg?branch=master)](https://travis-ci.org/snakemake-workflows/Phage_clouds_snake)

This is the snakemake implementation of the Phage Clouds workflow. At the moment, the workflow's structure is designed to work in the CBD server, however it can be customized for a more general use. It takes a directory containing genbank records of phage genome sequences, and returns an intergenomic distance network of phages whose clusters are defined based on user-defined host genus and distance threshold. The following diagram illustrates the overall Phage_clouds workflow:

![Phage_clouds](https://github.com/guille0387/Phage_clouds_snake/blob/master/dag.png)

Database icon taken from [www.flaticon.com](https://www.flaticon.com/)

## Authors

* Guillermo Rangel (@guille0387)

## Usage

### Step 1: Install Snakemake

It is recommended to install snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and creating an independent environment for it. Instructions on how to install snakemake can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

However, from within the CBD server it is possible to activate a ready-to-use environment that has snakemake installed in it. Simply run the following command from the terminal:

`conda activate ~/projects/snakemake`

### Step 2: Clone git repo and add input data

[Clone](https://help.github.com/en/articles/cloning-a-repository) this github repository to a location within your allocated
area in the CBD filesystem where you would like to conduct your analysis. This workflow was designed to take as input a directory containing all the phage genomes' genbank files that are to be part of the intergenomic distance network. Genbank files must have a .gbk extension, however these may be gzipped (.gz). It is recommended to either add the input directory to the root of the cloned repository, or to create a symlink to the input directory from within the repository's root.

This is a private repository, so if you don't have access to it please request access to Guillermo Rangel before attempting to clone it into the CBD server.

### Step 3: Configure workflow

In addition to taking as input the directory containing the genbank files, the workflow also requires that the user indicates a bacterial host genus they are interested in and a distance threshold to be applied to the mash distance values. All these parameters must be set using the file config.yaml, located within the config directory. Edit the file to add the following details:

input_dir: *name of directory containing gbk files*  
host: *bacterial host genus of your interest*  
threshold: *mash distance threshold to be applied to intergenomic distance network*

### Step 4: Execute workflow

Run the command below from within the repository's root to execute the workflow. **NOTE**: make sure you've activated the conda environment that contains snakemake.

`snakemake --cores 10 --use-conda`

### Step 5: Check results

All the output files generated during workflow execution will be saved in the results directory. Furthermore, the first time you run the workflow a copy of NCBI's taxonomy database will be downloaded and used to generate a .sqlite instance of it, both of which are kept in the resources directory.

Providing execution of the complete workflow was successful, the requested host-centric graph should be found in the results directory as *Host*_*threshold*_clouds.html. This file can be visualised in any wen browser such as Chrome, Safari, Firefox.
