#------------------------------------------------------------------------
# Title: Set Network Parameters for Various Cutoffs
# Date: May 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# set network parameters
if(percentage_cutoff == 0.2){
  # adipose parameters
  adipose_soft_power <- 4
  adipose_deep_split <- 3
  adipose_ME_cor_thres <- 0.1
  
  # liver parameters
  liver_soft_power <- 3
  liver_deep_split <- 3
  liver_ME_cor_thres <- 0.1
  
} else if(percentage_cutoff == 0.1){
  # adipose parameters
  adipose_soft_power <- 4
  adipose_deep_split <- 3
  adipose_ME_cor_thres <- 0.1
  
  # liver parameters
  liver_soft_power <- 4
  liver_deep_split <- 3
  liver_ME_cor_thres <- 0.1
} else{
  # adipose parameters
  adipose_soft_power <- 4
  adipose_deep_split <- 3
  adipose_ME_cor_thres <- 0.1
  
  # liver parameters
  liver_soft_power <- 4
  liver_deep_split <- 3
  liver_ME_cor_thres <- 0.1
}
