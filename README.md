# HAPPI_GWAS
A genome-wide association study (GWAS) tool written in R

## Usage
1. Clone the repository: git clone https://github.com/angelovicilab/HAPPI_GWAS.git
2. cd HAPPI_GWAS
3. Rscript setup.R
4. Rscript HAPPI_GWAS.R \<yaml file\> [COMMANDS]
   - COMMANDS:
     - -help (yaml file is not required)
     - -generateBLUP
     - -generateBLUE
     - -GAPIT
     - -extractHaplotype (optional, dependent on GAPIT parameter)
     - -searchGenes (optional, dependent on GAPIT parameter)

## Upgrade Packages
update.packages(ask=FALSE, checkBuilt = TRUE)

## Install Packages
<ul style="list-style-type:square">
  <li>install.packages("dplyr", dependencies = TRUE)</li>
  <li>install.packages("tidyr", dependencies = TRUE)</li>
  <li>install.packages("ggplot2", dependencies = TRUE)</li>
  <li>install.packages("tibble", dependencies = TRUE)</li>
  <li>install.packages("stringr", dependencies = TRUE)</li>
  <li>install.packages("gplots", dependencies = TRUE)</li>
  <li>install.packages("ape", dependencies = TRUE)</li>
  <li>install.packages("BiocManager", dependencies = TRUE)</li>
  <li>install.packages("car", dependencies = TRUE)</li>
  <li>install.packages("data.table", dependencies = TRUE)</li>
  <li>install.packages("DataCombine", dependencies = TRUE)</li>
  <li>install.packages("EMMREML", dependencies = TRUE)</li>
  <li>install.packages("foreach", dependencies = TRUE)</li>
  <li>install.packages("doParallel", dependencies = TRUE)</li>
  <li>install.packages("lme4", dependencies = TRUE)</li>
  <li>install.packages("scatterplot3d", dependencies = TRUE)</li>
  <li>install.packages("genetics", dependencies = TRUE)</li>
  <li>install.packages("LDheatmap", dependencies = TRUE)</li>
  <li>install.packages("gridExtra", dependencies = TRUE)</li>
  <li>install.packages("yaml", dependencies = TRUE)</li>
  <li>install.packages("bigmemory", dependencies = TRUE)</li>
  <li>install.packages("biganalytics", dependencies = TRUE)</li>
  <li>BiocManager::install("Biobase", update = TRUE, ask = FALSE)</li>
  <li>BiocManager::install("BiocGenerics", update = TRUE, ask = FALSE)</li>
  <li>BiocManager::install("snpStats", update = TRUE, ask = FALSE)</li>
  <li>BiocManager::install("multtest", update = TRUE, ask = FALSE)</li>
  <li>BiocManager::install("zlibbioc", update = TRUE, ask = FALSE)</li>
</ul>

## Updates
May 15, 2020: 
We have made changes to the tool and hosted it on https://github.com/Angelovici-Lab/HAPPI.GWAS. Please checkout the latest tool.
