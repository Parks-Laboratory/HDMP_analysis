#------------------------------------------------------------------------
# Title: Run WGCNA for HMDP Expression Data (Adipose and Liver)
# Date: May 17 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------


# CONTENTS:
  # Pull data from database and filter out genes
  # Run QC for missing data & outliers
  # Choose the soft-thresholding power - set params in separate wcgna_network_params.R
  # Run network (functions)
  # Format cluster data into data.frame
  # Save data

# NOTES:
  # filter: keep gene if at least 20% of strains show at least 20% up/down change in expression from chow to highfat
  # genes left over: adipose, liver
  # final data format: strains long, probe id wide
#------------------------------------------------------------------------


# Input Command Line Args
#------------------------------------------------------------------------

library(optparse)

# set arguments
option_list <- list(
  make_option(c("--sql"), default = FALSE, help = "Whether to run SQL to obtain data (only need to run once) [default %default]"),
  make_option(c("--percentage_cutoff"), default = 0.2, help = "Cutoff for up/down change in expression [default %default]"), 
  make_option(c("--output_folder"), help = "Folder for code outputs"), 
  make_option(c("--opt_softpower"), default = FALSE, help = "Whether to run softpower analysis [default %default]")
)

# process arguments
opt <- parse_args(OptionParser(option_list = option_list))


# Load Libraries
#------------------------------------------------------------------------
library(reshape2)
library(stringr)
library(plyr)
library(dplyr)
library(dtplyr)
library(tidyr)
library(sqldf)
library(data.table)
library(RODBC)
library(WGCNA)

# set input parameters
sql <- opt$sql
percentage_cutoff <- opt$percentage_cutoff
run_softpower <- opt$opt_softpower

setwd("E:/MICROARRAY DATA/REACTIVE_HDMP_NETWORKS/")
source("E:/MICROARRAY DATA/wgcna_pipeline/wgcna_network_pipeline.R")
source("E:/MICROARRAY DATA/REACTIVE_HDMP_NETWORKS/wcgna_network_params.R")

code_output_folder <- opt$output_folder
code_output_folder <- file.path(code_output_folder, percentage_cutoff)
if(!dir.exists(code_output_folder)) dir.create(code_output_folder)


# Process SQL Data
#------------------------------------------------------------------------

if(sql){
  
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
      subset(str_detect(TABLE_NAME, "expression_diff"))
    
    # generate the SQL queries from the table names
    queries <- paste("select * from", paste(data_sources$TABLE_SCHEM, data_sources$TABLE_NAME, sep = "."), "order by probe_id, strain")
    
    # connect to SQL Server and extract data
    data <- lapply(queries, function(q) sqlQuery(db, q))
    names(data) <- str_extract(queries, "adipose|liver")
    
    # return sql results
    return(data)
  }
  
  # function to find the percentage of strains that have percent change of gene expression 20% up or down by gene
  # Params: tissue, adipose or liver
  #         cutoff, percentage cutoff for up/down percent
  # Returns: data.frame of the percent changes
  grab_filter <- function(tissue, cutoff){
    
    # obtain the names of table/views from SQL Server
    # select only the tables that contain differences to run network algorithm on
    data_sources <- sqlTables(db) %>% 
      subset(TABLE_SCHEM == "dbo") %>% 
      dplyr::select(TABLE_CAT, TABLE_SCHEM, TABLE_NAME) %>% 
      subset(str_detect(TABLE_NAME, tissue) & str_detect(TABLE_NAME, "expression_by_strain") & str_detect(TABLE_NAME, "_male"))
    
    # generate the SQL queries from the table names
    queries <- paste("select * from", paste(data_sources$TABLE_SCHEM, data_sources$TABLE_NAME, sep = "."))
    
    # connect to SQL Server and extract data
    highfat <- sqlQuery(db, str_subset(queries, "highfat"))
    chow <- sqlQuery(db, str_subset(queries, "chow"))
    
    # run SQL query to find the percentage of strains that have percent change of gene expression either 20% up or down by gene
    results <- sqldf(paste0("
                            select probe_id, AVG(expression_percent_change > ", cutoff, " OR expression_percent_change < -", cutoff, ") as filter_option
                            from(
                            select highfat.gene_symbol, highfat.chr, highfat.probe_id, highfat.strain, (highfat.mean_rma_expression - chow.mean_rma_expression) / chow.mean_rma_expression as expression_percent_change
                            from highfat
                            inner join 
                            chow
                            on highfat.strain = chow.strain and highfat.probe_id = chow.probe_id
                            ) as tmp 
                            group by probe_id
                            order by filter_option DESC
                            "))
  
  # return the results of percent change
  return(results)
}

  # extract data
  diff_data <- grab_sql()
  
  # extract the filtering data
  filter_adipose <- grab_filter("adipose", percentage_cutoff)
  filter_liver <- grab_filter("liver", percentage_cutoff)
  
  # close the database
  odbcClose(db)
  
  
  # Filter Out Genes: keep if at least 20% of strains show at least 20% up/down change in expression from chow to highfat
  #------------------------------------------------------------------------
  
  # find the genes in which at least 20% of strains had big change in gene expression
  sub_filter_adipose <- subset(filter_adipose, filter_option >= 0.2)
  sub_filter_liver <- subset(filter_liver, filter_option >= 0.2)
  
  # QC: how many genes are left over after filter
  cat("how many genes are left over after filter:")
  nrow(sub_filter_adipose) 
  nrow(sub_filter_liver) 
  
  # filter difference data to keep the genes specified
  filter_data <- list()
  filter_data$adipose <- subset(diff_data$adipose, probe_id %in% sub_filter_adipose$probe_id)
  filter_data$liver <- subset(diff_data$liver, probe_id %in% sub_filter_liver$probe_id)
  
  
  # Process SQL Data
  #------------------------------------------------------------------------
  
  # process SQL data
  process_sql_data(sql_data = filter_data, use_col = "expression_difference", output_folder = code_output_folder)

}


# Load Data
#------------------------------------------------------------------------

load( file.path(code_output_folder, "network_expr_data.Rdata") ) # called expr_data and strain_names


# Choose the soft-thresholding power
#------------------------------------------------------------------------

if(run_softpower){
  pick_soft_power(expr_data[["adipose"]], output_folder = code_output_folder, output_suffix = "_adipose")
  pick_soft_power(expr_data[["liver"]], output_folder = code_output_folder, output_suffix = "_liver")
}


# Run network (code)
#------------------------------------------------------------------------

# import the affy annotation data
db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
affy <- sqlQuery(db, "select probesetID as probeid, gene_symbol, entrez_gene from dbo.affy_microarray_annotation")
odbcClose(db)

# import clinical trait differences
db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
clinical_traits <- sqlQuery(db, "select * from clinical_trait_differences_male")
odbcClose(db)

run_wgcna(
  data_id = "adipose", 
  expr_data = expr_data$adipose, 
  probeid = colnames(expr_data$adipose), 
  strains = strain_names$adipose, 
  affy = affy, 
  clinical_traits = clinical_traits,
  soft_power = adipose_soft_power, 
  min_mod_size = 30, 
  deep_split = adipose_deep_split, 
  ME_cor_thres = adipose_ME_cor_thres, 
  output_folder = code_output_folder,
  output_genetree = FALSE
)


run_wgcna(
  data_id = "liver", 
  expr_data = expr_data$liver, 
  probeid = colnames(expr_data$liver), 
  strains = strain_names$liver, 
  affy = affy,
  clinical_traits = clinical_traits,
  soft_power = liver_soft_power, 
  min_mod_size = 30, 
  deep_split = liver_deep_split, 
  ME_cor_thres = liver_ME_cor_thres, 
  output_folder = code_output_folder,
  output_genetree = FALSE
)
