# Metagenomics of Cheese Rind Microbiomes :cheese: :dna:

This repository contains the CSV files used for making the figures associated with the [Metagenomics of Cheese Rind Microbiomes pub](insert doi once it is available) as well as a Jupyter notebook with the code used to make the figures. Data from this pub, including raw reads and assemblies, can be found in the [European Nucleotide Archive under Project PRJEB58160](https://www.ebi.ac.uk/ena/browser/view/PRJEB58160). This includes a total of 24 metagenomic assemblies from five distinct washed rind cheese microbiomes. For three cheeses, paired short and long read sequencing data from three different timepoints in the cheese aging process is available.

The "Viral prediction summary tables" folder within 'processed_data' contains summaries of the combined viral prediction tables from [GeNomad 1.2.0 with database v1.1](https://zenodo.org/record/7015982#.Y9hoMi2B28c) and [CheckV 1.0.1 with database v1.4](https://www.nature.com/articles/s41587-020-00774-7) for nine of the long read metagenomic assemblies.

The "cheesegenomes-k31-scaled1k" within 'databases' folder contains sourmash signatures from previously sequenced microbial genomes from cheese. These signatures were used in the 'sourmash gather' analysis.
