# HAPPI_GWAS
A genome-wide association study (GWAS) tool written in R.

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
  <li>install.packages("dplyr")</li>
  <li>install.packages("tidyr")</li>
  <li>install.packages("ggplot2")</li>
  <li>install.packages("gplots")</li>
  <li>install.packages("ape")</li>
  <li>install.packages("BiocManager")</li>
  <li>install.packages("car")</li>
  <li>install.packages("data.table")</li>
  <li>install.packages("DataCombine")</li>
  <li>install.packages("EMMREML")</li>
  <li>install.packages("foreach")</li>
  <li>install.packages("doParallel")</li>
  <li>install.packages("lme4")</li>
  <li>install.packages("scatterplot3d")</li>
  <li>install.packages("genetics")</li>
  <li>install.packages("LDheatmap")</li>
  <li>install.packages("gridExtra")</li>
  <li>install.packages("yaml")</li>
</ul>

## Reminder
Please run git pull on terminal before open or run the HAPPI_GWAS pipeline since this application is still being developed.


