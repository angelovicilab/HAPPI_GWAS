## Clean the workspace
rm(list = ls())



# Set repository so that required packages can be downloaded
r = getOption("repos")
r["CRAN"] = "https://cran.cnr.berkeley.edu/"
options(repos = r)



R_library = paste(R.version$platform, "-library", sep = "")
R_version = gsub("\\.0$", "", paste(R.version$major, R.version$minor, sep = "."))

# Install path
p <- file.path("~/R", R_library, R_version)
if(!dir.exists(p)){
  dir.create(path = p, recursive = TRUE)
}



# Gather required packages
packages <- c("ape",
              "BiocManager",
              "compiler", "car",
              "data.table", "DataCombine",
              "EMMREML",
              "genetics", "gplots", "gridExtra",
              "lme4", "LDheatmap",
              "scatterplot3d",
              "dplyr", "tidyr", "tibble", "ggplot2", "stringr", 
              "yaml",
              "foreach", "doParallel")

# Check packages and install them if needed
invisible(lapply(packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "https://cran.cnr.berkeley.edu/", lib = p)
    library(x, lib.loc = p, character.only = TRUE)
  }
}))



# The packages here are from BiocManager
bioc_packages <- c("zlibbioc", "snpStats", "multtest")

# Check packages and install them if needed
invisible(lapply(bioc_packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, lib.loc = p, lib = p)
    library(x, lib.loc = p, character.only = TRUE)
  }
}))



# Print this after all packages are successfully installed
loaded_packages <- sessionInfo()
print(names(loaded_packages$otherPkgs))



# Import EMMA
source("http://www.zzlab.net/GAPIT/emma.txt")

# Import FarmCPU
source("http://www.zzlab.net/FarmCPU/FarmCPU_functions.txt")

# Import GAPIT
# source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
source("http://www.zzlab.net/GAPIT/previous/gapit_functions20191105.txt")

# Source R files
source("func_print_help.R")
source("func_read_file.R")
source("func_remove_duplicates.R")
source("func_outlier_removal.R")
source("func_boxcox_transformation.R")
source("func_generate_BLUP.R")
source("func_generate_BLUE.R")
source("func_farming_with_GAPIT.R")
source("func_extract_haplotype.R")
source("func_search_genes.R")



## Get command line arguments
args <- commandArgs(trailingOnly = TRUE)


#######################################################################
## Help support and extra configuration
#######################################################################

# help is in args
if (any(c("-h", "-help", "--help") %in% args)) {
  print_help()
  quit(status = 0)
}

# set number of cores
cores = 1
if ("-cores" %in% args) {
  index <- match("-cores", args)

  cores = ifelse(!is.null(index) & !is.na(index) & !is.na(as.integer(index)) & !is.na(as.integer(args[as.integer(index)+1])), as.integer(args[as.integer(index)+1]),1)
  cores = ifelse(!is.na(detectCores(logical = FALSE)) & cores > detectCores(logical = FALSE), detectCores(logical = FALSE), cores)

}



#######################################################################
## Starting and initialization
#######################################################################

cat(rep("\n", 2));print("-------------------- HAPPI_GWAS Start --------------------");cat(rep("\n", 2))

registerDoParallel(cores = cores)

cat(rep("\n", 2))
print(paste("Number of cores will be used is ", getDoParWorkers(), sep = ""))
cat(rep("\n", 2))



#######################################################################
## Check provided arguments and read in YAML file data
#######################################################################

# Check if the YAML file path is in the args
if (!identical(args, character(0)) & length(args) > 0 & file.exists(file.path(args[1]))) {
  print("YAML file exists!!!")
  cat(rep("\n", 2))
} else{
  print("The first provided argument is not a YAML file path or the YAML file does not exists!!!")
  cat(rep("\n", 2))
  print_help()
  quit(status = -1)
}

# If there are more than one argument, print all the actions needed to perform
if (length(args) > 1) {
  print(args[2:length(args)])
  cat(rep("\n", 2))
} else{
  print("No parameter! No action is required to perform!!!")
  cat(rep("\n", 2))
  print_help()
  quit(status = -1)
}


## Import configuration file and read in raw data and reference files
yaml_dat <- tryCatch({
  read_yaml(file.path(args[1]))
}, error = function(e) {
  print("The yaml file is invalid!!!")
  cat(rep("\n", 2))
  print_help()
  quit(status = -1)
})



#######################################################################
## Check yaml data parameters
#######################################################################

if (exists("yaml_dat")) {
  if(is.null(yaml_dat$raw_data) & is.null(yaml_dat$BLUP)){
    print("The yaml data does not have input file path!!!")
    quit(status = -1)
  } 
  if(is.null(yaml_dat$by_column) | is.null(yaml_dat$start_column) | 
      is.null(yaml_dat$BLUP_by_column) | is.null(yaml_dat$BLUP_start_column)){
        print("The yaml data does not have input file by_column or start_column configuration!!!")
        quit(status = -1)
  }
  if(is.null(yaml_dat$output)){
    print("The yaml data does not have output path!!!")
    quit(status = -1)
  }
} else{
  print("The yaml data was not read into the work space!!!")
  quit(status = -1)
}



#######################################################################
## Read in all file and data that are specified in the YAML file
#######################################################################

cat(rep("\n", 2))
## Import raw data
raw_data <- read_file(file_path = yaml_dat$raw_data)
if (is.null(raw_data)) {
  print("The raw_data parameter is NULL.")
} else{
  print("raw_data has been loaded into memory.")
}

if (!is.null(yaml_dat$by_column)) {
  by_column <- yaml_dat$by_column
  print(paste("by_column: ", by_column, sep = ""))
} else{
  print("The by_column parameter is NULL.")
}

if (!is.null(yaml_dat$start_column)) {
  start_column <- yaml_dat$start_column
  print(paste("start_column: ", start_column, sep = ""))
} else{
  print("The start_column parameter is NULL.")
}

cat(rep("\n", 2))
## Import BULP data
BLUP <- read_file(file_path = yaml_dat$BLUP)
if (is.null(BLUP)) {
  print("The BLUP parameter is NULL.")
} else{
  print("BLUP has been loaded into memory.")
}

if (!is.null(yaml_dat$BLUP_by_column)) {
  BLUP_by_column <- yaml_dat$BLUP_by_column
  print(paste("BLUP_by_column: ", BLUP_by_column, sep = ""))
} else{
  print("The BLUP_by_column parameter is NULL.")
}

if (!is.null(yaml_dat$BLUP_start_column)) {
  BLUP_start_column <- yaml_dat$BLUP_start_column
  print(paste("BLUP_start_column: ", BLUP_start_column, sep = ""))
} else{
  print("The BLUP_start_column parameter is NULL.")
}

cat(rep("\n", 2))
## Import GAPIT reference files
# GAPIT kinship matrix
GAPIT_kinship_matrix <- read_file(file_path = yaml_dat$GAPIT_kinship_matrix, header = FALSE)
if (is.null(GAPIT_kinship_matrix)) {
  print("The GAPIT_kinship_matrix parameter is NULL.")
} else{
  print("GAPIT_kinship_matrix has been loaded into memory.")
}

# GAPIT covariates
GAPIT_covariates <- read_file(file_path = yaml_dat$GAPIT_covariates)
if (is.null(GAPIT_covariates)) {
  print("The GAPIT_covariates parameter is NULL.")
} else{
  print("GAPIT_covariates has been loaded into memory.")
}

# GAPIT hapmap file
GAPIT_hapmap <- read_file(file_path = yaml_dat$GAPIT_hapmap, header = FALSE)
if (is.null(GAPIT_hapmap)) {
  print("The GAPIT_hapmap parameter is NULL.")
} else{
  print("GAPIT_hapmap has been loaded into memory.")
}

# GAPIT genotype data (numeric)
GAPIT_genotype_data_numeric <- read_file(file_path = yaml_dat$GAPIT_genotype_data_numeric)
if (is.null(GAPIT_genotype_data_numeric)) {
  print("The GAPIT_genotype_data_numeric parameter is NULL.")
} else{
  print("GAPIT_genotype_data_numeric has been loaded into memory.")
}

# GAPIT genotype map (numeric)
GAPIT_genotype_map_numeric <- read_file(file_path = yaml_dat$GAPIT_genotype_map_numeric)
if (is.null(GAPIT_genotype_map_numeric)) {
  print("The GAPIT_genotype_map_numeric parameter is NULL.")
} else{
  print("GAPIT_genotype_map_numeric has been loaded into memory.")
}

# GAPIT hapmap file extension
if (!is.null(yaml_dat$GAPIT_hapmap_file_extension)) {
  GAPIT_hapmap_file_extension <- yaml_dat$GAPIT_hapmap_file_extension
  print(paste("GAPIT_hapmap_file_extension: ", GAPIT_hapmap_file_extension, sep = ""))
} else{
  GAPIT_hapmap_file_extension <- NULL
  print("The GAPIT_hapmap_file_extension parameter is NULL.")
}

# GAPIT genotype data numeric file extension
if (!is.null(yaml_dat$GAPIT_genotype_data_numeric_file_extension)) {
  GAPIT_genotype_data_numeric_file_extension <- yaml_dat$GAPIT_genotype_data_numeric_file_extension
  print(paste("GAPIT_genotype_data_numeric_file_extension: ", GAPIT_genotype_data_numeric_file_extension, sep = ""))
} else{
  GAPIT_genotype_data_numeric_file_extension <- NULL
  print("The GAPIT_genotype_data_numeric_file_extension parameter is NULL.")
}

# GAPIT genotype map numeric file extension
if (!is.null(yaml_dat$GAPIT_genotype_map_numeric_file_extension)) {
  GAPIT_genotype_map_numeric_file_extension <- yaml_dat$GAPIT_genotype_map_numeric_file_extension
  print(paste("GAPIT_genotype_map_numeric_file_extension: ", GAPIT_genotype_map_numeric_file_extension, sep = ""))
} else{
  GAPIT_genotype_map_numeric_file_extension <- NULL
  print("The GAPIT_genotype_map_numeric_file_extension parameter is NULL.")
}

# GAPIT hapmap filename
if (!is.null(yaml_dat$GAPIT_hapmap_filename)) {
  GAPIT_hapmap_filename <- yaml_dat$GAPIT_hapmap_filename
  print(paste("GAPIT_hapmap_filename: ", GAPIT_hapmap_filename, sep = ""))
} else{
  GAPIT_hapmap_filename <- NULL
  print("The GAPIT_hapmap_filename parameter is NULL.")
}

# GAPIT genotype data numeric filename
if (!is.null(yaml_dat$GAPIT_genotype_data_numeric_filename)) {
  GAPIT_genotype_data_numeric_filename <- yaml_dat$GAPIT_genotype_data_numeric_filename
  print(paste("GAPIT_genotype_data_numeric_filename: ", GAPIT_genotype_data_numeric_filename, sep = ""))
} else{
  GAPIT_genotype_data_numeric_filename <- NULL
  print("The GAPIT_genotype_data_numeric_filename parameter is NULL.")
}

# GAPIT genotype map numeric filename
if (!is.null(yaml_dat$GAPIT_genotype_map_numeric_filename)) {
  GAPIT_genotype_map_numeric_filename <- yaml_dat$GAPIT_genotype_map_numeric_filename
  print(paste("GAPIT_genotype_map_numeric_filename: ", GAPIT_genotype_map_numeric_filename, sep = ""))
} else{
  GAPIT_genotype_map_numeric_filename <- NULL
  print("The GAPIT_genotype_map_numeric_filename parameter is NULL.")
}

# GAPIT genotype file path
if (!is.null(yaml_dat$GAPIT_genotype_file_path)) {
  GAPIT_genotype_file_path <- file.path(yaml_dat$GAPIT_genotype_file_path)
  if(dir.exists(file.path(GAPIT_genotype_file_path))){
    print(paste("GAPIT_genotype_file_path: ", GAPIT_genotype_file_path, sep = ""))
  } else{
    print("The GAPIT_genotype_file_path parameter is a file path that does not exists.")
    quit(status = -1)
  }
} else{
  GAPIT_genotype_file_path <- NULL
  print("The GAPIT_genotype_file_path parameter is NULL.")
}

# GAPIT genotype file named sequentially from
if (!is.null(yaml_dat$GAPIT_genotype_file_named_sequentially_from)) {
  GAPIT_genotype_file_named_sequentially_from <- as.numeric(yaml_dat$GAPIT_genotype_file_named_sequentially_from)
  print(paste("GAPIT_genotype_file_named_sequentially_from: ", GAPIT_genotype_file_named_sequentially_from, sep = ""))
} else{
  GAPIT_genotype_file_named_sequentially_from <- 0
  print("The GAPIT_genotype_file_named_sequentially_from parameter is 0.")
}

# GAPIT genotype file named sequentially to
if (!is.null(yaml_dat$GAPIT_genotype_file_named_sequentially_to)) {
  GAPIT_genotype_file_named_sequentially_to <- as.numeric(yaml_dat$GAPIT_genotype_file_named_sequentially_to)
  print(paste("GAPIT_genotype_file_named_sequentially_to: ", GAPIT_genotype_file_named_sequentially_to, sep = ""))
} else{
  GAPIT_genotype_file_named_sequentially_to <- 0
  print("The GAPIT_genotype_file_named_sequentially_to parameter is 0.")
}

# GAPIT model
if (!is.null(yaml_dat$GAPIT_model)) {
  GAPIT_model <- yaml_dat$GAPIT_model
  print(paste("GAPIT_model: ", GAPIT_model, sep = ""))
} else{
  GAPIT_model <- NULL
  print("The GAPIT_model parameter is NULL.")
}

# GAPIT SNP.MAF
if (!is.null(yaml_dat$GAPIT_SNP_MAF)) {
  GAPIT_SNP_MAF <- as.numeric(yaml_dat$GAPIT_SNP_MAF)
  if(is.na(GAPIT_SNP_MAF)) { GAPIT_SNP_MAF <- 0 }
  print(paste("GAPIT_SNP_MAF: ", GAPIT_SNP_MAF, sep = ""))
} else{
  GAPIT_SNP_MAF <- 0
  print("The GAPIT_SNP_MAF parameter is 0.")
}

# GAPIT PCA total
if (!is.null(yaml_dat$GAPIT_PCA_total)) {
  GAPIT_PCA_total <- as.numeric(yaml_dat$GAPIT_PCA_total)
  if(is.na(GAPIT_PCA_total)) { GAPIT_PCA_total <- 0 }
  print(paste("GAPIT_PCA_total: ", GAPIT_PCA_total, sep = ""))
} else{
  GAPIT_PCA_total <- 0
  print("The GAPIT_PCA_total parameter is 0.")
}

# GAPIT Model selection
if (!is.null(yaml_dat$GAPIT_Model_selection)) {
  GAPIT_Model_selection <- as.logical(yaml_dat$GAPIT_Model_selection)
  if(!is.logical(GAPIT_Model_selection) | is.na(GAPIT_Model_selection)) { GAPIT_Model_selection <- FALSE }
  print(paste("GAPIT_Model_selection: ", GAPIT_Model_selection, sep = ""))
} else{
  GAPIT_Model_selection <- FALSE
  print("The GAPIT_Model_selection parameter is FALSE.")
}

# GAPIT SNP test
if (!is.null(yaml_dat$GAPIT_SNP_test)) {
  GAPIT_SNP_test <- as.logical(yaml_dat$GAPIT_SNP_test)
  if(!is.logical(GAPIT_SNP_test) | is.na(GAPIT_SNP_test)) { GAPIT_SNP_test <- FALSE }
  print(paste("GAPIT_SNP_test: ", GAPIT_SNP_test, sep = ""))
} else{
  GAPIT_SNP_test <- FALSE
  print("The GAPIT_SNP_test parameter is FALSE.")
}

# GAPIT file output
if (!is.null(yaml_dat$GAPIT_file_output)) {
  GAPIT_file_output <- as.logical(yaml_dat$GAPIT_file_output)
  if(!is.logical(GAPIT_file_output) | is.na(GAPIT_file_output)) { GAPIT_file_output <- FALSE }
  print(paste("GAPIT_file_output: ", GAPIT_file_output, sep = ""))
} else{
  GAPIT_file_output <- FALSE
  print("The GAPIT_file_output parameter is FALSE.")
}

# GAPIT p.Value threshold
if (!is.null(yaml_dat$GAPIT_p_value_threshold)) {
  GAPIT_p_value_threshold <- as.numeric(yaml_dat$GAPIT_p_value_threshold)
  if(is.logical(GAPIT_p_value_threshold)){ GAPIT_p_value_threshold <- NA }
  print(paste("GAPIT_p_value_threshold: ", GAPIT_p_value_threshold, sep = ""))
} else{
  GAPIT_p_value_threshold <- NA
  print("The GAPIT_p_value_threshold parameter is NA.")
}

# GAPIT p.Value FDR threshold
if (!is.null(yaml_dat$GAPIT_p_value_fdr_threshold)) {
  GAPIT_p_value_fdr_threshold <- as.numeric(yaml_dat$GAPIT_p_value_fdr_threshold)
  if(is.logical(GAPIT_p_value_fdr_threshold)){ GAPIT_p_value_fdr_threshold <- NA }
  print(paste("GAPIT_p_value_fdr_threshold: ", GAPIT_p_value_fdr_threshold, sep = ""))
} else{
  GAPIT_p_value_fdr_threshold <- NA
  print("The GAPIT_p_value_fdr_threshold parameter is NA.")
}

# GAPIT LD_number
if (!is.null(yaml_dat$GAPIT_LD_number)) {
  GAPIT_LD_number <- as.numeric(yaml_dat$GAPIT_LD_number)
  print(paste("GAPIT_LD_number: ", GAPIT_LD_number, sep = ""))
} else{
  print("The GAPIT_LD_number parameter is NULL.")
}

cat(rep("\n", 2))
## Haploview
# Haploview_file_path
if (!is.null(yaml_dat$Haploview_file_path)) {
  Haploview_file_path <- normalizePath(file.path(yaml_dat$Haploview_file_path))
  if(dir.exists(file.path(Haploview_file_path))){
    print(paste("Haploview_file_path: ", Haploview_file_path, sep = ""))
  } else{
    print("The Haploview_file_path parameter is a file path that does not exists.")
    quit(status = -1)
  }
} else{
  Haploview_file_path <- NULL
  print("The Haploview_file_path parameter is NULL.")
}

# Haploview_file_name
if (!is.null(yaml_dat$Haploview_file_name)) {
  Haploview_file_name <- as.character(yaml_dat$Haploview_file_name)
  print(paste("Haploview_file_name: ", Haploview_file_name, sep = ""))
} else{
  Haploview_file_name <- NULL
  print("The Haploview_file_name parameter is NULL.")
}

# Haploview_file_extension
if (!is.null(yaml_dat$Haploview_file_extension)) {
  Haploview_file_extension <- as.character(yaml_dat$Haploview_file_extension)
  print(paste("Haploview_file_extension: ", Haploview_file_extension, sep = ""))
} else{
  Haploview_file_extension <- NULL
  print("The Haploview_file_extension parameter is NULL.")
}

# Haploview_file_named_sequentially_from
if (!is.null(yaml_dat$Haploview_file_named_sequentially_from)) {
  Haploview_file_named_sequentially_from <- as.numeric(yaml_dat$Haploview_file_named_sequentially_from)
  print(paste("Haploview_file_named_sequentially_from: ", Haploview_file_named_sequentially_from, sep = ""))
} else{
  Haploview_file_named_sequentially_from <- NA
  print("The Haploview_file_named_sequentially_from parameter is NA.")
}

# Haploview_file_named_sequentially_to
if (!is.null(yaml_dat$Haploview_file_named_sequentially_to)) {
  Haploview_file_named_sequentially_to <- as.numeric(yaml_dat$Haploview_file_named_sequentially_to)
  print(paste("Haploview_file_named_sequentially_to: ", Haploview_file_named_sequentially_to, sep = ""))
} else{
  Haploview_file_named_sequentially_to <- NA
  print("The Haploview_file_named_sequentially_to parameter is NA.")
}

cat(rep("\n", 2))
## Match Gene Start and Gene Stop
# GFF_file_path
if (!is.null(yaml_dat$GFF_file_path)) {
  GFF_file_path <- normalizePath(file.path(yaml_dat$GFF_file_path))
  if(dir.exists(file.path(GFF_file_path))){
    print(paste("GFF_file_path: ", GFF_file_path, sep = ""))
  } else{
    print("The GFF_file_path parameter is a file path that does not exists.")
    quit(status = -1)
  }
} else{
  GFF_file_path <- NULL
  print("The GFF_file_path parameter is NULL.")
}

# GFF_file_name
if (!is.null(yaml_dat$GFF_file_name)) {
  GFF_file_name <- as.character(yaml_dat$GFF_file_name)
  print(paste("GFF_file_name: ", GFF_file_name, sep = ""))
} else{
  GFF_file_name <- NULL
  print("The GFF_file_name parameter is NULL.")
}

# GFF_file_extension
if (!is.null(yaml_dat$GFF_file_extension)) {
  GFF_file_extension <- as.character(yaml_dat$GFF_file_extension)
  print(paste("GFF_file_extension: ", GFF_file_extension, sep = ""))
} else{
  GFF_file_extension <- NULL
  print("The GFF_file_extension parameter is NULL.")
}

# GFF_file_named_sequentially_from
if (!is.null(yaml_dat$GFF_file_named_sequentially_from)) {
  GFF_file_named_sequentially_from <- as.numeric(yaml_dat$GFF_file_named_sequentially_from)
  print(paste("GFF_file_named_sequentially_from: ", GFF_file_named_sequentially_from, sep = ""))
} else{
  GFF_file_named_sequentially_from <- NA
  print("The GFF_file_named_sequentially_from parameter is NA.")
}

# GFF_file_named_sequentially_to
if (!is.null(yaml_dat$GFF_file_named_sequentially_to)) {
  GFF_file_named_sequentially_to <- as.numeric(yaml_dat$GFF_file_named_sequentially_to)
  print(paste("GFF_file_named_sequentially_to: ", GFF_file_named_sequentially_to, sep = ""))
} else{
  GFF_file_named_sequentially_to <- NA
  print("The GFF_file_named_sequentially_to parameter is NA.")
}

cat(rep("\n", 2))
## Create output folder
if (!is.null(yaml_dat$output)) {
  if (!dir.exists(file.path(yaml_dat$output))) {
    dir.create(path = file.path(yaml_dat$output), showWarnings = TRUE, recursive = TRUE)
    if (dir.exists(file.path(yaml_dat$output))) {
      output <- normalizePath(file.path(yaml_dat$output))
      print(paste0("The output folder is created. Output path: ", output))
    } else{
      print("The output folder cannot be created.")
      quit(status = -1)
    }
  } else{
    output <- normalizePath(file.path(yaml_dat$output))
    print(paste0("The output folder exists. Output path: ", output))
  }
} else{
  print("The output parameter is NULL.")
  quit(status = -1)
}

cat(rep("\n", 2))

#######################################################################
## Run action base on arguments
#######################################################################

# # removeDuplicates
# if (all("-removeDuplicates" %in% args)) {
#   index <- match("-removeDuplicates", args)
#   print(paste(index, ": removeDuplicates", sep = ""))
#
#   if (exists("raw_data") & exists("by_column") & exists("start_column") & dir.exists(output)) {
#
#     folder_path <- file.path(output, "removeDuplicates")
#
#     if (!dir.exists(folder_path)) {
#       dir.create(path = folder_path, showWarnings = TRUE, recursive = TRUE)
#     } else{
#       print("The removeDuplicates folder exists.")
#     }
#
#     # Using customized function to remove duplicates
#     raw_data <- remove_duplicates(dat = raw_data, by_column = by_column)
#
#     write.csv(x = raw_data, file = file.path(folder_path, "raw_data.csv"), row.names = TRUE, na = "")
#   }
# }
#
# # outlierRemoval
# if (all("-outlierRemoval" %in% args)) {
#   index <- match("-outlierRemoval", args)
#   print(paste(index, ": outlierRemoval", sep = ""))
#
#   if (exists("raw_data") & !is.null(raw_data) & exists("by_column") & exists("start_column") & dir.exists(output)) {
#
#     folder_path <- file.path(output, "outlierRemoval")
#
#     if (!dir.exists(folder_path)) {
#       dir.create(path = folder_path, showWarnings = TRUE, recursive = TRUE)
#     } else{
#       print("The outlierRemoval folder exists.")
#     }
#
#     # Using customized function to remove outliers
#     results <- outlier_removal(dat = raw_data, by_column = by_column, start_column = start_column)
#
#     if (is.list(results)) {
#       raw_data <- results$Outlier_removed_data
#
#       capture.output( results$Outliers_residuals, file = file.path(folder_path, "Outliers_residuals.txt"))
#       write.csv(x = results$Outlier_data, file = file.path(folder_path, "Outlier_data.csv"), row.names = FALSE, na = "" )
#       write.csv( x = results$Outlier_removed_data, file = file.path(folder_path, "Outlier_removed_data.csv"), row.names = FALSE, na = "")
#     } else{
#       raw_data <- results
#
#       write.csv(x = results, file = file.path(folder_path, "Outlier_removed_data.csv"), row.names = FALSE, na = "")
#       print("No outlier has been found!!!")
#     }
#
#   }
# }
#
# # boxcoxTransformation
# if (all("-boxcoxTransformation" %in% args)) {
#   index <- match("-boxcoxTransformation", args)
#   print(paste(index, ": boxcoxTransformation", sep = ""))
#
#   if (exists("raw_data") & !is.null(raw_data) & exists("by_column") & exists("start_column") & dir.exists(output)) {
#
#     folder_path <- file.path(output, "boxcoxTransformation")
#
#     if (!dir.exists(folder_path)) {
#       dir.create(path = folder_path, showWarnings = TRUE, recursive = TRUE)
#     } else{
#       print("The boxcoxTransformation folder exists.")
#     }
#
#     # Using customized function to perform box-cox transformation
#     results <- boxcox_transformation(dat = raw_data, by_column = by_column, start_column = start_column)
#
#     if (is.list(results)) {
#       raw_data <- results$Boxcox_transformed_data
#
#       write.csv(x = results$Lambda_values, file = file.path(folder_path, "Lambda_values.csv"), row.names = FALSE, na = "")
#       write.csv(x = results$Boxcox_transformed_data, file = file.path(folder_path, "Boxcox_transformed_data.csv"), row.names = FALSE, na = "")
#     } else{
#       raw_data <- results
#
#       print("No lambda found!!! Data returned without transformed!!!")
#     }
#
#   }
# }

# generateBLUP
if (all("-generateBLUP" %in% args)) {
  index <- match("-generateBLUP", args)
  print(paste(index, ": generateBLUP", sep = ""))

  if (exists("raw_data") & !is.null(raw_data) & exists("by_column") & exists("start_column") & dir.exists(output)) {

    folder_path <- file.path(output, "generateBLUP")

    if (!dir.exists(folder_path)) {
      dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
    } else{
      print("The generateBLUP folder exists.")
    }

    # Using customized function to generate BLUP
    results <-
      generate_BLUP(dat = raw_data, by_column = by_column, start_column = start_column)

    if (is.list(results)) {
      BLUP <- results$BLUP
      BLUP_by_column <- 1
      BLUP_start_column <- 2

      write.csv( x = results$BLUP, file = file.path(folder_path, "BLUP.csv"), row.names = FALSE, na = "" )
      write.csv(x = results$Lambda_values, file = file.path(folder_path, "Lambda_values.csv"), row.names = FALSE, na = "")
      write.csv(x = results$Boxcox_transformed_data, file = file.path(folder_path, "Boxcox_transformed_data.csv"), row.names = FALSE, na = "")
      capture.output( results$Outliers_residuals, file = file.path(folder_path, "Outliers_residuals.txt"))
      write.csv(x = results$Outlier_data, file = file.path(folder_path, "Outlier_data.csv"), row.names = FALSE, na = "" )
      write.csv( x = results$Outlier_removed_data, file = file.path(folder_path, "Outlier_removed_data.csv"), row.names = FALSE, na = "")
    } else{
      print("No BLUP generated!!!")
    }

  }
}

# generateBLUE
if (all("-generateBLUE" %in% args)) {
  index <- match("-generateBLUE", args)
  print(paste(index, ": generateBLUE", sep = ""))

  if (exists("raw_data") & !is.null(raw_data) & exists("by_column") & exists("start_column") & dir.exists(output)) {

    folder_path <- file.path(output, "generateBLUE")

    if (!dir.exists(folder_path)) {
      dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
    } else{
      print("The generateBLUE folder exists.")
    }

    # Using customized function to generate BLUP
    results <-
      generate_BLUE(dat = raw_data, by_column = by_column, start_column = start_column)

    if (is.list(results)) {
      BLUP <- results$BLUE
      BLUP_by_column <- 1
      BLUP_start_column <- 2

      write.csv( x = results$BLUE, file = file.path(folder_path, "BLUE.csv"), row.names = FALSE, na = "" )
      write.csv(x = results$Lambda_values, file = file.path(folder_path, "Lambda_values.csv"), row.names = FALSE, na = "")
      write.csv(x = results$Boxcox_transformed_data, file = file.path(folder_path, "Boxcox_transformed_data.csv"), row.names = FALSE, na = "")
      capture.output( results$Outliers_residuals, file = file.path(folder_path, "Outliers_residuals.txt"))
      write.csv(x = results$Outlier_data, file = file.path(folder_path, "Outlier_data.csv"), row.names = FALSE, na = "" )
      write.csv( x = results$Outlier_removed_data, file = file.path(folder_path, "Outlier_removed_data.csv"), row.names = FALSE, na = "")
    } else{
      print("No BLUE generated!!!")
    }

  }
}

# GAPIT
if (all("-GAPIT" %in% args)) {
  index <- match("-GAPIT", args)
  print(paste(index, ": GAPIT", sep = ""))

  if (exists("BLUP") & !is.null(BLUP) & exists("BLUP_by_column") & exists("BLUP_start_column") &
      exists("GAPIT_LD_number") & GAPIT_LD_number >= 0 & dir.exists(output)) {

    folder_path <- file.path(output, "GAPIT")

    if (!dir.exists(folder_path)) {
      dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
    } else{
      print("The GAPIT folder exists.")
    }

    # Using customized function to run GAPIT
    combined_gwas_result <-
      farming_with_GAPIT(
        dat = BLUP,
        by_column = BLUP_by_column,
        start_column = BLUP_start_column,
        output_path = folder_path,
        p_value_threshold = GAPIT_p_value_threshold,
        p_value_fdr_threshold = GAPIT_p_value_fdr_threshold,
        ld_number = GAPIT_LD_number,
        KI = GAPIT_kinship_matrix,
        CV = GAPIT_covariates,
        G = GAPIT_hapmap,
        GD = GAPIT_genotype_data_numeric,
        GM = GAPIT_genotype_map_numeric,
        file.Ext.G = GAPIT_hapmap_file_extension,
        file.Ext.GD = GAPIT_genotype_data_numeric_file_extension,
        file.Ext.GM = GAPIT_genotype_map_numeric_file_extension,
        file.G = GAPIT_hapmap_filename,
        file.GD = GAPIT_genotype_data_numeric_filename,
        file.GM = GAPIT_genotype_map_numeric_filename,
        file.path = GAPIT_genotype_file_path,
        file.from = GAPIT_genotype_file_named_sequentially_from,
        file.to = GAPIT_genotype_file_named_sequentially_to,
        model = GAPIT_model,
        SNP.MAF = GAPIT_SNP_MAF,
        PCA.total = GAPIT_PCA_total,
        Model.selection = GAPIT_Model_selection,
        SNP.test = GAPIT_SNP_test,
        file.output = GAPIT_file_output
      )

  }
}

# Extract Haplotype
if (all("-extractHaplotype" %in% args)) {
  index <- match("-extractHaplotype", args)
  print(paste(index, ": extractHaplotype", sep = ""))

      # if combined_gwas_result is not null and other requirements are satisfied then extract haplotype
      if(exists("combined_gwas_result") & !is.null(combined_gwas_result) & !is.null(Haploview_file_path) & !is.null(Haploview_file_name) &
          !is.null(Haploview_file_extension) & !is.na(Haploview_file_named_sequentially_from) & !is.na(Haploview_file_named_sequentially_to)){

            folder_path <- file.path(output, "GAPIT")

            if (!dir.exists(folder_path)) {
              dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
            } else{
              print("The GAPIT folder exists.")
            }

            combined_gwas_result <- extract_haplotype(
              combined_gwas_result = combined_gwas_result,
              output_path = folder_path,
              Haploview_file_path = Haploview_file_path,
              Haploview_file_name = Haploview_file_name,
              Haploview_file_extension = Haploview_file_extension,
              Haploview_file_named_sequentially_from = Haploview_file_named_sequentially_from,
              Haploview_file_named_sequentially_to = Haploview_file_named_sequentially_to
            )
      }
}

# searchGenes
if (all("-searchGenes" %in% args)) {
  index <- match("-searchGenes", args)
  print(paste(index, ": searchGenes", sep = ""))

      # if combined_gwas_result is not null and other requirements are satisfied then search genes
      if(exists("combined_gwas_result") & !is.null(combined_gwas_result) & !is.null(GFF_file_path) & !is.null(GFF_file_name) &
          !is.null(GFF_file_extension) & !is.na(GFF_file_named_sequentially_from) & !is.na(GFF_file_named_sequentially_to)){

            folder_path <- file.path(output, "GAPIT")

            if (!dir.exists(folder_path)) {
              dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
            } else{
              print("The GAPIT folder exists.")
            }

            combined_gwas_result <- search_genes(
              combined_gwas_result = combined_gwas_result,
              output_path = folder_path,
              GFF_file_path = GFF_file_path,
              GFF_file_name = GFF_file_name,
              GFF_file_extension = GFF_file_extension,
              GFF_file_named_sequentially_from = GFF_file_named_sequentially_from,
              GFF_file_named_sequentially_to = GFF_file_named_sequentially_to
            )
      }
}



#######################################################################
## Ending and releasing resources
#######################################################################

cat(rep("\n", 2));print("-------------------- HAPPI_GWAS Exit --------------------");cat(rep("\n", 2))

stopImplicitCluster()

