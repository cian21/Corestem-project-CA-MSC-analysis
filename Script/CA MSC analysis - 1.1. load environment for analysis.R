# Corestem project - CA MSC Analysis.proj
# 20201215 Sehwan Chun at Corestem, Inc.
# 1.1. load environment for analysis

#### 1. Library Loading ####
packs = c("DESeq2","org.Hs.eg.db","data.table","ggpubr", "edgeR")
lapply(packs, require, character.only = TRUE)
rm(packs)

#### 2. Files Loading ####
expr_file = "./Data/Rawdata/gene_count_matrix.csv"
expr_file = read.csv(expr_file, stringsAsFactors = F, header = T)
expr_file = expr_file[complete.cases(expr_file),]

group_file = "./Data/Rawdata/ca_group.txt"
group_file = read.table(group_file, stringsAsFactors = F, header = T)

save.image(file = "./Data/CA MSC analysis files.image")
