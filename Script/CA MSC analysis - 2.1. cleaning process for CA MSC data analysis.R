# Corestem project - CA MSC Analysis.proj
# 20201215 Sehwan Chun at Corestem, Inc.
# 2.1. cleaning process for CA MSC data analysis

#### 1. source Loading ####
load("./Data/CA MSC analysis files.image")
source("./Script/Functions/CA MSC analysis - functions.R")

#### 2. cpm cleaning ####
row.names(expr_file) = expr_file$gene_id
expr_file$gene_id = NULL

expr_cpm = CPM_extract_run(expr_file,group_file)
expr_DEG = DEG_extract_run(expr_file,group_file)

tmp = group_file[c(5:9,1:4,10:13),3] #for order
diff_list = NULL
diff_nums = NULL
for (i in 1:nrow(expr_cpm)){
    if(cor.test(as.numeric(expr_cpm[i,1:13]),tmp)$p.value < 0.05){
        diff_list = c(diff_list, as.character(expr_cpm[i,15]))
        diff_nums = c(diff_nums, as.character(expr_cpm[i,14]))
    }
}

expr_DEG_filtered = subset(expr_DEG, Symbol %in% setdiff(expr_DEG$Symbol,diff_list))
expr_DEG_filtered$FDR = p.adjust(expr_DEG_filtered$PValue, method = "fdr")

save.image(file = "./Data/CA MSC analysis files cleaned.image")
