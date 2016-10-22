#------------------------------------------------------------------------
# Title: Process CEL files for HMDP expression data
# Date: April 8 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS:
  # Download/Open GEO Files
  # Normalization with RMA
  # Open Annotation Files 
  # Combine Annotation Files with ExpressionSet Data
  # Upload to SQL server

# NOTES:
  # run the download of CEL files and normalization on HTCondor (more RAM) 
  # save normalized results as Rdata for easier loading
#------------------------------------------------------------------------


# Load Libraries
#------------------------------------------------------------------------
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(oligo)
library(RODBC)


# Download/Open GEO Files
#------------------------------------------------------------------------
# set the GEO file name
GEO_file <- c("GSE64770", "GSE42890", "GSE16780")

# open microarry data
setwd("E:/MICROARRAY DATA/REACTIVE_HDMP_NETWORKS/")

# download GEO files and unpack it - only need to do this once
# lapply(GEO_file, GEOquery::getGEOSuppFiles)
# lapply(GEO_file, function(x) untar(paste0(x, "/", x, "_RAW.tar"), exdir = "CEL"))

# extract names of celfiles
celfiles <- list.files("CEL", full.names = TRUE)
sub_celfiles <- str_subset(celfiles, "(cel|CEL).gz")

# open up cel files
# rawData <- read.celfiles(sub_celfiles)


# Normalization using RMA
#------------------------------------------------------------------------
# perform rma
# normData <- rma(rawData)


# Save Results
#------------------------------------------------------------------------
# process takes a long time - just save results as Rdata files
# save(normData, file = "normData.Rdata")

# load the data
load("E:/MICROARRAY DATA/HDMP_raw_cel_files/normData.Rdata")


# Undo Log2 Transformation (in the rma)
#------------------------------------------------------------------------
# exprs(normData) <- 2^exprs(normData)


# Open Annotation Files from GEO and Lab
#------------------------------------------------------------------------

# open file containing translated highfat annotations from BP
experimental_data_highfat <- fread("Annotation/ORIGINAL_HMDP_HIGHFAT_IDS_STRAINS_ANNOTATION.csv")
setnames(experimental_data_highfat, tolower(colnames(experimental_data_highfat)))

# convert the mouse id into just the id string (to match the annotations from GEO)
setnames(experimental_data_highfat, "mouse id", "mouse_id")
experimental_data_highfat <- dplyr::mutate(experimental_data_highfat, mouse_id = str_extract(mouse_id, "\\d+"))

# open file containing translated chow annotations from BP
experimental_data_chow <- fread("Annotation/ORIGINAL_HMDP_CHOW_IDS_STRAINS_ANNOTATION.csv")
setnames(experimental_data_chow, tolower(colnames(experimental_data_chow)))
setnames(experimental_data_chow, c("mouseid", "jackson_name"), c("mouse_id", "strain"))

# open files containing annotations from GEO into a list of files
file_paths <- list.files("GEO_accession", full.names = TRUE)
files <- lapply(file_paths, function(x){
  
  # read in file, select needed cols, add file name
  f <- fread(x)
  setnames(f, tolower(colnames(f)))
  f <- dplyr::select(f, accession, title)
  f <- dplyr::mutate(f, source = str_replace_all(x, "GEO_accession/|_accession.csv", ""))
  setnames(f, "title", "mouse_id")
  
  # return opened file
  return(f)
})
names(files) <- str_extract(file_paths, "adipose|liver|HMDP")

# manually edit (column names, string formatting) 3 types of files so the match each other
files[["adipose"]] <- tidyr::separate(files[["adipose"]], source, c("diet", "tissue"), sep = "_")
files[["liver"]] <- tidyr::separate(files[["liver"]], source, c("diet", "tissue"), sep = "_")
files[["liver"]] <- dplyr::mutate(files[["liver"]], mouse_id = str_replace(mouse_id, "UCLA_mouse_liver_", ""))
files[["HMDP"]] <- files[["HMDP"]] %>% 
  tidyr::separate(mouse_id, c("mouse_id", "tissue"), sep = "_") %>% 
  dplyr::mutate(diet = "highfat", source = NULL)


# Merge Annotation Files
#------------------------------------------------------------------------

# CHOW ADIPOSE: merge annotations 
files[["adipose"]] <- merge(files[["adipose"]], experimental_data_chow, "mouse_id", all.x = TRUE) %>% 
  dplyr::select(-number)

# obtain the strain information from chow data 
files[["liver"]] <- files[["liver"]] %>% 
  dplyr::mutate(
    number = str_extract(mouse_id, "\\d+$"),
    mouse_id = str_extract(mouse_id, "[a-zA-Z\\d-]*_\\d+")
  ) 

# recode a mouse id so that it is mergeable with the experimental data annotations
files[["liver"]] <- dplyr::mutate(files[["liver"]], mouse_id = ifelse(mouse_id == "B6Cc3-1_6", "B6Cc3_6", mouse_id))

# CHOW LIVER: merge annotations 
files[["liver"]] <- merge(files[["liver"]], experimental_data_chow, "mouse_id", all.x = TRUE) %>% 
  dplyr::select(-matches("number"))

# HIGHFAT: merge annotations 
files[["HMDP"]] <- merge(files[["HMDP"]], experimental_data_highfat, "mouse_id", all.x = TRUE)

# reformat mouse id to bet the same width
files[["HMDP"]] <- mutate(files[["HMDP"]], mouse_id = str_pad(mouse_id, width = 4, side = "left", pad = "0"))

# combine 3 files
GEO_annotate <- data.table::rbindlist(files, use.names = TRUE, fill = TRUE)


# Annotate ExpressionSet
#------------------------------------------------------------------------
# edit the names of normData to match the Accession (the first sequence of characters)
sampleNames <- str_extract(sampleNames(normData), "[:alnum:]+")
pData(normData)$accession <- sampleNames

# add the phenotype data into the expression set
pData(normData) <- merge(pData(normData), GEO_annotate, "accession", all.x = TRUE) 
sampleNames(normData) <- sampleNames


# Format for Database
#------------------------------------------------------------------------

# connect to DB
db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')

# function to upload specified diet and tissue into database as its own table
upload_to_db <- function(in_diet, in_tissue){
  
  # find the indices of the samples to keep coresponding to diet and tissue
  keep <- with(pData(normData), diet == in_diet & tissue == in_tissue)
  sub_data <- subset(pData(normData), diet == in_diet & tissue == in_tissue)
  
  # subset and keep the samples that correspond to parameters
  subNormData <- normData[, keep]
  
  # melt the data from the expression set
  db_ready <- melt( exprs(subNormData) )
  colnames(db_ready) <- c("probe_id", "mouse_id", "rma_expression")

  # replace the GEO accession with the mouse id
  db_ready <- dplyr::mutate(db_ready, mouse_id = plyr::mapvalues(mouse_id, sub_data$accession, sub_data$mouse_id))

  # upload data to the database
  sqlSave(channel = db, dat = db_ready, tablename = paste(in_diet, in_tissue, "microarray", sep = "_"), rownames = FALSE)  
  
  # return table that was uploaded to DB
  return(db_ready)
}

# run uploads into sql database
Sys.time()
x1 <- upload_to_db("chow", "liver")
x2 <- upload_to_db("chow", "adipose")
x3 <- upload_to_db("highfat", "liver")
x4 <- upload_to_db("highfat", "adipose")
Sys.time() # 6 hours!

odbcClose(db)
