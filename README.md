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
</ul>

## Reminder
Please run git pull on terminal before open or run the HAPPI_GWAS pipeline since this application is still being developed.


