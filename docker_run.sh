#!/bin/bash

yaml_file=$1
gapit=$2
extract_haplotype=$3
search_genes=$4

cd /home/HAPPI_GWAS/ && Rscript HAPPI_GWAS.R $yaml_file $gapit $extract_haplotype $search_genes
