#------------------------------------------------------------------------
# Title: Format Network Data for Cytoscape
# Date: May 31 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS:
  # Pull data (expression data, network data, sql data - for color)
  # Subset to Modules
  # Format for Cytoscape

# NOTES:
#------------------------------------------------------------------------



# Input Command Line Args
#------------------------------------------------------------------------

library(optparse)

# set arguments
option_list <- list(
  make_option(c("--tissue"), default = "liver", help = "What tissue?"), 
  make_option(c("--module"), help = "What module color?"), 
  make_option(c("--threshold"), default = 0.9, help = "Threshold for inclusion [0-1]")
)

# process arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set input parameters
modules <- opt$module
threshold <- opt$threshold
tissue <- opt$tissue

# input/output
setwd("E:/MICROARRAY DATA/REACTIVE_HDMP_NETWORKS/code_outputs/0.2_with_0.1_threshold/")
output_folder <- paste0("cytoscape/thres", threshold, "/")

if(!dir.exists(output_folder)) dir.create(output_folder)

# Load Libraries
#------------------------------------------------------------------------
library(reshape2)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(RODBC)
library(sqldf)
library(WGCNA)

source("E:/MICROARRAY DATA/wgcna_pipeline/wgcna_network_pipeline.R")


# Pull Expression Data
#------------------------------------------------------------------------

load(paste0("E:/MICROARRAY DATA/REACTIVE_HDMP_NETWORKS/code_outputs/0.2_with_0.1_threshold/network_expr_data.Rdata")) # called expr_data
expr_data <- expr_data[[tissue]]


# Pull Network Data
#------------------------------------------------------------------------

load(paste0(tissue, "_networkConstruction.Rdata"))


# Pull Data from Database
#------------------------------------------------------------------------

# connect to SQL Server
db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')

# function to find the direction of expression change from chow to highfat
# Params: tissue, adipose or liver
# Returns: data.frame of the percent changes
grab_filter <- function(tissue){
  
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
  
  # run SQL query to find the direction of expression change from chow to highfat (up or down)
  results <- sqldf(paste0("          
                          select probe_id, gene_symbol, SUM(direction_change > 0) - SUM(direction_change < 0) as majority_direction_change
                          from(
                            select highfat.gene_symbol, highfat.chr, highfat.probe_id, sign(highfat.mean_rma_expression - chow.mean_rma_expression) as direction_change
                            from highfat
                            inner join 
                            chow
                            on highfat.strain = chow.strain and highfat.probe_id = chow.probe_id
                          ) as tmp 
                          group by probe_id
                          "))
  
  # return the results of percent change
  return(results)
}

# find the direction of expression change from chow to highfat (up or down)
dir_change <- grab_filter(tissue = tissue)

# close the database
odbcClose(db)


# Subset to Interesting modules
#------------------------------------------------------------------------

# extract dir change to: probeid, gene_symbo, direction change majority
dir_change <- dir_change %>% 
  mutate(probeid = probe_id, change = ifelse(majority_direction_change > 0, "up", "down")) %>% 
  dplyr::select(probeid, change) %>% 
  distinct()

# combine direction change and network data
network_data$networkDF <- merge(network_data$networkDF, dir_change, "probeid", all = TRUE)


# Format for Cytoscape
#------------------------------------------------------------------------

# export the network into edge and node list files Cytoscape can read
format_cytoscape(network_data = network_data, expr_data = expr_data, modules = modules, threshold = threshold, include_col = "change", output_folder = output_folder, output_suffix = tissue)

# import files to edit
edges <- fread( paste0(output_folder, "CytoscapeInput-edges-", paste(modules, collapse=" - "), ".txt") )
nodes <- fread( paste0(output_folder, "CytoscapeInput-nodes-", paste(modules, collapse=" - "), ".txt") )
nodes <- mutate(nodes, fromNode = nodeName, nodeName = NULL, altName = NULL, change = direction, direction = NULL)

# combine node/edge files
network <- merge(edges, nodes, "fromNode")

# export network file for Cytoscope
fwrite(network, paste0(output_folder, "Cytoscape_", modules, ".txt"))
