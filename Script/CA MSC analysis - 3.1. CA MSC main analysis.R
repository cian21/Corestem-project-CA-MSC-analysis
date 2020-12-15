# Corestem project - CA MSC Analysis.proj
# 20201215 Sehwan Chun at Corestem, Inc.
# 3.1. CA MSC main analysis

#### 1. source Loading ####
load("./Data/CA MSC analysis files cleaned.image")
source("./Script/Functions/CA MSC analysis - functions.R")

#### 2. confirm DEGs and GO ####
expr_DEG_filtered = subset(expr_DEG_filtered, FDR < 0.1)

#just for GO analysis, the dispersion was changed by exclusion of genes 
expr_DEG_filtered_goana = DEG_extract_filtered_run(expr_file, group_file, diff_nums) 