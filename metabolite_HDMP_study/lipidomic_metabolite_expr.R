#' Compute Differential Lipid Metabolite Presence for Control vs. Agpat5 Mice

rm(list= ls())

# Load Libraries 
#----------------------------------------------------------------------------------------------------
library(readxl)
library(reshape2)
library(stringr)
library(plyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(gtools)
library(dplyr)
theme_set(theme_bw())


# Open Raw Data 
#----------------------------------------------------------------------------------------------------
# read in raw data in excel format
raw_lipids_pos <- read_excel("C:/Users/jnnguyen2/Downloads/PosMode_Data_053015_TripleTof.xlsx")
raw_lipids_neg <- read_excel("C:/Users/jnnguyen2/Downloads/Copy_of_Parks_NegMode_071215_Results.xlsx")


# Pre-Processing 
#--------------------------------------------------------------------------------------------------
clean_raw_data <- function(raw_data){
  
  # select raw data columns (positions 2-27) 
  sub_data <- data.table(raw_data[, c(2:27)])
  
  # rename sub_data columns 
  colnames(sub_data)[1] <- "lipids_samples_periods" # original name: Lipids\\Samples-Periods
  colnames(sub_data)[2] <- "experimental_details_count" # original name: Experiment Details: Count
  
  # unite the two columns from above to obtain a unique experiment/metabolite name 
  sub_data <- unite(sub_data, metabolite, lipids_samples_periods, experimental_details_count, sep = " ")
  
  # replace the old column names (formatted "ParksPos##-#, period:1") with shorter names ("Pos##-#")
  old_names <- str_subset(colnames(sub_data), "Parks(Pos|Neg)\\d+")
  new_names <- str_extract(old_names, "(Pos|Neg)\\d+-\\d+")
  setnames(sub_data, old_names, new_names)
  
  # melt sub_data to obtain the mouse # and duplicate number long
  sub_data_expr <- melt(sub_data, variable.name = "mouse_id", value.name = "expression")
  return(sub_data_expr)
  
}

# clean data
clean_lipids_pos <- clean_raw_data(raw_lipids_pos)
clean_lipids_neg <- clean_raw_data(raw_lipids_neg)

# combine cleaned data
clean_lipids <- rbindlist(list(clean_lipids_pos, clean_lipids_neg))

# edit mouse id information
clean_lipids <- mutate(clean_lipids, mouse_id = str_extract(mouse_id, "\\d+-\\d+"))
clean_lipids <- mutate(clean_lipids, mouse_id = ifelse(nchar(mouse_id) == 3, paste0("0", mouse_id), mouse_id))


# Crosswalk Mouse to Treatment
#----------------------------------------------------------------------------------------------------

# generate side-by-side the mouse_id and the treatment it received
mouse_num <- c("02", "07", "08", "10", "04", "05", "06", "09", "01", "03", "11", "12")
trt_grp <- c(rep("AGPAT5-ASO", 4), rep("Control-ASO", 4), rep("Saline", 4))
xwalk <- data.frame(mouse_num, trt_grp, stringsAsFactors = FALSE)


# Function to Process Data
# Inputs:
  # clean lipids: data frames with columns metabolite, mouse id (with duplicate info), expressions
  # xwalk: a crosswalk of mouse id to treatment
  # rm_trt: remove any treatment ("Saline")
  # opt_rm_8: specific to application, removing mouse 8 with low knockout expression
#----------------------------------------------------------------------------------------------------

process_data <- function(clean_lipids, xwalk, rm_trt = "Saline", opt_rm_8 = FALSE){
  
  
  # Add Grouping Identifiers 
  #--------------------------------------------------------------------------------------------------
  
  # generate a mouse id label (# 1-12) from the rownames and removing the duplication number
  lipid_data_all <- separate(clean_lipids, mouse_id, c("mouse_id", "dup"))

  # generate a new variable with the treatment
  lipid_data_all <- mutate(lipid_data_all, trt = mapvalues(mouse_id, xwalk$mouse_num, xwalk$trt_grp))
  
  
  # QC: Assessing Variability in Duplicates
  #--------------------------------------------------------------------------------------------------
  
  # set the duplicates wide in side by side columns 
  qc_duplicates <- lipid_data_all %>% 
    dcast(metabolite + mouse_id + trt ~ dup, value.var = "expression")
  
  # add the prefix "dup" before the duplicate number
  setnames(qc_duplicates, as.character(1:2), paste0("dup", 1:2))
  
  # compute the means of the duplicates
  row_means <- qc_duplicates %>% 
    dplyr::select(matches("dup")) %>% 
    apply(1, mean)
  
  # compute and order by the absolute difference between two duplicates
  qc_duplicates <- mutate(qc_duplicates, diff = abs(dup1 - dup2), mean = row_means) 
  
  # find the metabolite in which the duplicates varied greatly for trt & control
  flag_metabolites <- qc_duplicates %>% 
    # remove saline
    subset(trt != "Saline") %>% 
    # find the metabolite whose duplicate differences were off by more than 10 for at least 4 mice
    group_by(metabolite) %>% 
    dplyr::summarise(n_diff = sum(diff > 100)) %>% 
    subset(n_diff >= 4)
  flag_metabolites <- flag_metabolites$metabolite
  
  
  # QC: remove the 8th mouse frm the AGPAT group as a potential outlier
  #--------------------------------------------------------------------------------------------------
  
  if(opt_rm_8) lipid_data_all <- subset(lipid_data_all, !str_detect(mouse_id, "8"))
  
  
  # Fix the 0 Issue
  #--------------------------------------------------------------------------------------------------
  # find the smallest expression that is not 0
  sorted_expressions <- sort(unique(lipid_data_all$expression))
  
  # divide that small expression level by 10 to add to other values
  min_val_to_add <- sorted_expressions[2]
  
  # add small value to everything
  lipid_data_expr_shift <- mutate(lipid_data_all, shift_expr = 0.001 + expression)
  

  # QC: Fix the Dubious Intensity Values (low signal to noise ratio)
  #--------------------------------------------------------------------------------------------------
  
  # remove saline
  lipid_data_expr_shift <- subset(lipid_data_expr_shift, trt != rm_trt)
  
  # across all tretatment and mice, compute summary statistics for the expression
  qc_intensity <- lipid_data_expr_shift %>% 
    group_by(metabolite) %>% 
    dplyr::summarise(
      mn = mean(expression), 
      max = max(expression), 
      min = min(expression), 
      iqr = IQR(expression), 
      range_diff = abs(diff(range(expression)))
    )
  
  # Potential Criteria: max < 10

  # remove those metabolites whose intensity is less than 10 (these are ones with measurement errors)
  omit <- subset(qc_intensity, max < 10)$metabolite

  # generate a message to tell us how much we lost
  b4 <- nrow(lipid_data_expr_shift)
  lipid_data_expr_shift <- subset(lipid_data_expr_shift, !(metabolite %in% omit))
  after <- nrow(lipid_data_expr_shift)

  notice <- paste0("We lose ", b4 - after, " data points.", " (which is about ", (b4 - after)/after, ")")
  message(notice)
  
  
  # Average Duplicates 
  #--------------------------------------------------------------------------------------------------
  
  # grouping by the mouse, compute the mean expression level for each gene
  lipids_no_dups <- lipid_data_expr_shift %>% 
    group_by(metabolite, trt, mouse_id) %>% 
    dplyr::summarise(avr_dup_expr = mean(shift_expr))
 
  
  # Obtain P-Values 
  #--------------------------------------------------------------------------------------------------
  
  # obtain the two-sample t.test p-value for AGPAT vs Control
  t_test_p <- lipids_no_dups %>%
    group_by(metabolite) %>% 
    # error control for t-test in which all values are the same (metabolite abundance is 0 across all trt)
    dplyr::summarise(
      p.val = try_default(expr = t.test(avr_dup_expr ~ trt)$p.val, default = 1, quiet = TRUE)
    )

                
  # Average Across Trt Groups 
  #--------------------------------------------------------------------------------------------------
  
  # grouping by treatment, compute the mean expression level for each gene
  lipids_by_trt <- lipids_no_dups %>% 
    group_by(trt, metabolite) %>% 
    dplyr::summarise(mn_expr_by_trt = mean(avr_dup_expr))
  

  # Examine the Fold Change of AGPAT5 over Control
  #--------------------------------------------------------------------------------------------------
  
  # cast experimental details wide
  lipids_by_trt <- dcast(lipids_by_trt, trt~metabolite, value.var = "mn_expr_by_trt")
  
  # order by treatment alphabetically
  lipids_by_trt <- arrange(lipids_by_trt, trt)
  
  # a function to compute the fold change of agpat expression to control expression
  f <- function(x) foldchange(x[1], x[2])
  
  # generate expression ratios for each metabolite (some formatting steps at the end)
  expr_ratios <- lipids_by_trt %>% 
    dplyr::select(-trt) %>% 
    summarise_each(funs(f)) %>% 
    melt(variable.name = "metabolite", value.name = "fold.change")
  
  
  # Combine the Information of P-values and Expression Ratios
  #--------------------------------------------------------------------------------------------------

  summarized_data <- merge(t_test_p, expr_ratios, "metabolite", all = TRUE)
  
  
  # Return Expression Ratios & QC
  #--------------------------------------------------------------------------------------------------
  
  r <- list(
      summarized_data = summarized_data, 
      cleaned_raw_data = lipid_data_all, 
      qc = list(qc_duplicates = qc_duplicates, qc_intensity = qc_intensity, 
                qc_min_value = min_val_to_add, qc_metabolites_dups = flag_metabolites)
    )
  
  return(r)
}


# Run Function on Data
#----------------------------------------------------------------------------------------------------

opt_8 <- function(clean_lipids, opt_rm_8 = TRUE){
  
  # suffix: with or without 8
  suffix <- ifelse(opt_rm_8, " without 8", " with 8")
  
  # run analysis pipeline on raw data
  processed_data <- process_data(clean_lipids, xwalk, opt_rm_8 = opt_rm_8, rm_trt = "Saline")

  # extract data
  all_expr_data <- processed_data$summarized_data
  all_clean_raw_data <- processed_data$cleaned_raw_data
  
  # flag metabolites with variable duplicates
  flag_metabolites <- processed_data$qc$qc_metabolites_dups
  all_expr_data <- mutate(all_expr_data, flag_metabolites = ifelse(metabolite %in% flag_metabolites, 1, 0))
  
  
  # Export CSVs
  #----------------------------------------------------------------------------------------------------
  
  write.csv(all_expr_data, paste0("C:/Users/jnnguyen2/Downloads/out/results", suffix, ".csv"), row.names = FALSE)
  write.csv(all_clean_raw_data, paste0("C:/Users/jnnguyen2/Downloads/out/raw_data", suffix, ".csv"), row.names = FALSE)
  
  
  # Plot Volcano Plot
  #----------------------------------------------------------------------------------------------------
  
  all_expr_data <- mutate(all_expr_data, flag_fc = ifelse(abs(fold.change) >= 2, "flag", "no flag"))
  
  g2 <- ggplot(data = all_expr_data, aes(foldchange2logratio(fold.change), -log10(p.val), color = flag_fc)) +
    geom_point(size = 1.5, alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), color = "darkred", linetype = 2) +
    scale_colour_manual(values = c("flag" = "darkred", "no flag" = "darkgrey")) +
    xlab("log2 fold change") +
    ylab("-log10 p.value") +
    ggtitle(paste0("P-Value vs. Fold-Change", suffix)) +
    ylim(-.1, 4.8) +
    theme(legend.position = "none")
  
  png(paste0("C:/Users/jnnguyen2/Downloads/out/Volcano_Plot", suffix, ".png"), height = 600, width = 600)
  print(g2)
  dev.off()
  
  
  # QC
  #----------------------------------------------------------------------------------------------------
  
  if(!opt_rm_8) {
    qc_duplicates <- processed_data$qc$qc_duplicates
    # write.csv(qc_duplicates, "C:/Users/jnnguyen2/Downloads/out/qc_results.csv", row.names = FALSE)
  }
  
  # return(all_expr_data)
}

# run analysis on AGPAT5 vs control with and without 8
wo8 <- opt_8(clean_lipids, opt_rm_8 = TRUE) 
w8 <- opt_8(clean_lipids, opt_rm_8 = FALSE) 

# run analysis on control ASO vs saline
control_vs_aso <- process_data(clean_lipids, xwalk, rm_trt = "AGPAT5-ASO")$summarized_data
write.csv(control_vs_aso, "C:/Users/jnnguyen2/Downloads/out/control_results.csv")





# Distribution Plots by Treatment
#----------------------------------------------------------------------------------------------------

raw_data <- fread("C:/Users/jnnguyen2/Downloads/out/raw_data with 8.csv")
results <- fread("C:/Users/jnnguyen2/Downloads/out/results without 8.csv")

# subset to metabolites with p-values less than 0.05
signif_metabolites <- subset(results, p.val <= 0.05)
other_metabolites <- subset(results, str_detect(metabolite, "LPS|LPI") | str_detect(metabolite, "^PA"))
other_metabolites <- arrange(other_metabolites, p.val)[1:6]
signif_metabolites <- rbind(signif_metabolites, other_metabolites)

sub_lipids <- merge(raw_data, signif_metabolites, "metabolite")

# remove the prefix in the mouse id
sub_lipids <- mutate(sub_lipids, mouse_id = str_extract(mouse_id, "\\d+"))

# clean up some formatting for plotting
sub_lipids <- sub_lipids %>% 
  mutate(
    mouse_id = ifelse(nchar(mouse_id) == 1, paste0("0", mouse_id), mouse_id), 
    x.var = paste0(trt, ": ", mouse_id)
  ) %>% 
  subset(trt != "Saline") %>% 
  subset(mouse_id != "08")

# rearrange order
sub_lipids <- arrange(sub_lipids, p.val, metabolite)

# extract loop metabolites
loop_metabolites <- unique(sub_lipids$metabolite)

pdf("C:/Users/jnnguyen2/Downloads/out/boxplots.pdf")

l_ply(loop_metabolites, function(m){
  
  # subset to specific metabolite
  plot_data <- subset(sub_lipids, metabolite == m) 

  plot_data$trt <- factor(plot_data$trt)
  plot_data$trt <- reorder(plot_data$trt, ifelse(plot_data$trt == "AGPAT5-ASO", 2, 1))
  
  # compute the mean expression
  # plot_data <- plot_data %>%
  #   group_by(x.var, trt) %>%
  #   summarise(expression = mean(expression), p.val = unique(p.val), fold.change = unique(fold.change), flag = unique(flag_metabolites))
  # 
  label <- paste0("p.val = ", round(unique(plot_data$p.val), 3), "\n", "fold.change = ", round(unique(plot_data$fold.change), 3))
  
  g <- ggplot(plot_data, aes(x = trt, y = expression)) +
    geom_boxplot()+
    # geom_bar(stat = "identity") +
    # geom_point(size = 2) +
    xlab("Treatment") + ylab("Intensity") + ggtitle(paste0(m, "\n", label)) + 
    theme(legend.position = "none")
    # theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  if(unique(plot_data$flag) == 1) print("This metabolite has highly variable duplicates")
  print(g)
  
})
dev.off()
