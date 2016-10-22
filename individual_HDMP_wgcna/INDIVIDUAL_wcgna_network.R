#------------------------------------------------------------------------
# Title: Run WGCNA for Individual Expression Data
# Date: May 23 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS:
  # Choose the soft-thresholding power
  # Run network (functions)
  # Format cluster data into data.frame
  # Save data

# NOTES:
  # files too large - need to split up and run on clusters
  # Files to be transferred to clusters
    # wgcna pipeline files
    # affy_data.csv
    # network_expr_data.Rdata
#------------------------------------------------------------------------


# Input Command Line Args
#------------------------------------------------------------------------

library(optparse)

# set arguments
option_list <- list(
  make_option(c("--sql"), default = FALSE, help = "Whether the code is extracting data from SQL [default %default]"),
  make_option(c("--condor"), default = FALSE, help = "Whether code is running on Condor [default %default]"), 
  make_option(c("--i"), help = "What iteration to run, ranges from 0-5"), 
  make_option(c("--output_folder"), help = "Where to save output"),
  make_option(c("--opt_softpower"), default = FALSE, help = "Whether to run softpower analysis [default %default]")
)

# process arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set input parameters
sql <- opt$sql
condor <- opt$condor
i <- as.numeric(opt$i) + 1
code_output_folder <- opt$output_folder
run_softpower <- opt$opt_softpower

sp <- 3

# Load Libraries
#------------------------------------------------------------------------
library(reshape2)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(WGCNA)


# Run Networks
#------------------------------------------------------------------------

# run on condor for more memory options
if(sql){

  library(sqldf)
  library(RODBC)
  
  # load network code
  source("E:/MICROARRAY DATA/wgcna_pipeline/wgcna_network_pipeline.R")
  
  # set working directory
  setwd("E:/MICROARRAY DATA/INDIVIDUAL_NETWORKS/")
  
  
  # Pull Data from Database
  #------------------------------------------------------------------------
  
  # connect to SQL Server
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
  
  # function to extract highfat v chow expression difference data
  # Returns: list of data from sql for different tissues
  grab_sql <- function(){
    
    # obtain the names of table/views from SQL Server: tables that contain differences to run networks
    data_sources <- sqlTables(db) %>% 
      subset(TABLE_SCHEM == "dbo") %>% 
      dplyr::select(TABLE_CAT, TABLE_SCHEM, TABLE_NAME) %>% 
      subset(str_detect(TABLE_NAME, "avg_expression_by_strain")) %>% 
      mutate(diet = str_extract(TABLE_NAME, "chow|highfat"), tissue = str_extract(TABLE_NAME, "adipose|liver"), gender = str_extract(TABLE_NAME, "female|male"))
    
    # generate the SQL queries from the table names
    queries <- paste("select * from", paste(data_sources$TABLE_SCHEM, data_sources$TABLE_NAME, sep = "."), "order by probe_id, strain")
    
    # connect to SQL Server and extract data
    data <- lapply(queries, function(q) sqlQuery(db, q))
    names(data) <- with(data_sources, paste(tissue, diet, gender, sep = "_"))
    
    # return sql results
    return(data)
  }
  
  # extract data
  data <- grab_sql()
  
  # close the database
  odbcClose(db)
  
  
  # Process SQL Data
  #------------------------------------------------------------------------
  
  # process SQL data
  process_sql_data(sql_data = filter_data, use_col = "mean_rma_expression", output_folder = "code_outputs")
  

  # Softpower Analysis
  #------------------------------------------------------------------------

  # load data
  load( file.path(code_output_folder, "network_expr_data.Rdata") ) # called expr_data and strain_names
  options <- names(expr_data)
  
  # pick_softpower
  if(run_softpower){
    
    for(opt in options){
      pick_soft_power(expr_data[[opt]], output_folder = code_output_folder, output_suffix = paste0("_", opt))
    }
    
  }

  
  # Obtain Affy Data
  #------------------------------------------------------------------------
  
  db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
  affy <- sqlQuery(db, "select probesetID as probeid, gene_symbol, entrez_gene from dbo.affy_microarray_annotation")
  odbcClose(db)
  
  write.csv(affy, file.path(code_output_folder, "affy_data.csv"), row.names = FALSE)
  
  
} else if(condor){
  
  # Load Data
  #------------------------------------------------------------------------
  
  # load network code
  source("wgcna_functions_condor.R")
  source("wgcna_network_pipeline.R")
  
  # load gene annotations
  affy <- fread("affy_data.csv")
  
  # load data
  kind <- options[i]
  load("network_expr_data.Rdata") # expr_data & strain_names
  options <- names(expr_data)
  
  # set params
  kind <- options[i]
  expr_data <- expr_data[[kind]]
  strain_names <- strain_names[[kind]]
  
  
  # Run Pipeline
  #------------------------------------------------------------------------
  
  run_wgcna(
    data_id = kind, 
    expr_data = expr_data, 
    probeid = colnames(expr_data), 
    strains = strain_names, 
    affy = affy, 
    # clinical_traits = , NO CLINICAL TRAITS
    
    soft_power = sp,
    min_mod_size = 30,
    deep_split = 3,
    ME_cor_thres = 0.1,
    
    output_folder = code_output_folder,
    output_genetree = TRUE
  )

  
# second portion of code requires db access - run on server
} else{
  
  # load network code
  source("E:/MICROARRAY DATA/wgcna_pipeline/wgcna_network_pipeline.R")
  
  # set working directory
  setwd("E:/MICROARRAY DATA/INDIVIDUAL_NETWORKS/")
  
  # load the data: expression data and Rdata from Condor
  load(file.path(code_output_folder, "network_expr_data.Rdata"))
  options <- names(expr_data)
  
  for(kind in options){
    
    # load data
    load(file.path(code_output_folder, paste0(kind, "_outdata.Rdata")))
    load(file.path(code_output_folder, paste0(kind, "_networkConstruction.Rdata")))
    # load(file.path(code_output_folder, paste0(kind, "_gene_module_data.Rdata")))
    
    # print merged modules
    png(file.path(code_output_folder, paste0("clustering_after_merge_", kind, ".png")), height = 500, width = 700)
    plotDendroAndColors(genes$geneTree, cbind(genes$dynamicColors, modules$mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    dev.off()
    
  }
  
}

