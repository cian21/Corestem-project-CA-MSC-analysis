# # Corestem project - CA MSC Analysis.proj
# 20201215 Sehwan Chun at Corestem, Inc.
# functions.R

#### 1. Library Loading ####
packs = c("DESeq2","org.Hs.eg.db","data.table","ggpubr", "edgeR")
lapply(packs, require, character.only = TRUE)
rm(packs)

#### 2. DEG analysis ####
anno_DEG = function(DGE_table){
    entrezid = row.names(DGE_table)
    cols = c("ENTREZID", "SYMBOL")
    geneList = biomaRt::select(org.Hs.eg.db,
                               keys = entrezid,
                               columns = cols,
                               keytype = "ENTREZID")
    DGE_table$entrizid = entrezid
    
    #File Cleaning
    for (i in 1:nrow(DGE_table)){
        if(is.na(geneList[i,2]) == F){
            row.names(DGE_table)[i] = geneList[i,2]
        }
    }
    DGE_table$Symbol = row.names(DGE_table)
    return(DGE_table)
}
CPM_extract_run = function(expr_file, group_file){
    
    DGE = DGEList(counts = expr_file)
    keep = filterByExpr(DGE)
    DGE = DGE[keep, , keep.lib.sizes = FALSE] 
    
    DGE = calcNormFactors(DGE)
    DGE = estimateDisp(DGE)
    
    DGECPM = cpm(DGE, log = FALSE)
    DGECPM = as.data.frame(DGECPM)
    DGECPM = anno_DEG(DGECPM)
    
    return(DGECPM)
}
DEG_extract_run = function(expr_file, group_file){
    expr_file = expr_file[,-c(6:9)]
    group_file = group_file[-c(1:4),]
    rownames(group_file) = 1:nrow(group_file)
    
    design_matrix = model.matrix(~Group, data = group_file)
    
    DGE = DGEList(counts = expr_file, group = group_file$Group)
    keep = filterByExpr(DGE)
    DGE = DGE[keep, , keep.lib.sizes = FALSE] 
    DGE = calcNormFactors(DGE)
    DGE = estimateDisp(DGE, design_matrix)
    
    DGETable = exactTest(DGE)
    DGETable = DGETable$table
    DGETable$FDR = p.adjust(DGETable$PValue, method =  "fdr")
    DGETable = anno_DEG(DGETable)
    return(DGETable)
}
DEG_extract_filtered_run = function(expr_file, group_file, exclude_list){
    expr_file = expr_file[,-c(6:9)]
    expr_file = subset(expr_file, row.names(expr_file) %in% setdiff(row.names(expr_file), exclude_list))
    group_file = group_file[-c(1:4),]
    rownames(group_file) = 1:nrow(group_file)
    
    design_matrix = model.matrix(~Group, data = group_file)
    
    DGE = DGEList(counts = expr_file, group = group_file$Group)
    keep = filterByExpr(DGE)
    DGE = DGE[keep, , keep.lib.sizes = FALSE] 
    DGE = calcNormFactors(DGE)
    DGE = estimateDisp(DGE, design_matrix)
    
    DGETable = exactTest(DGE)
    go = goana(DGETable)
    top_up_10 = topGO(go, ont="BP", sort="Up", n=10)
    top_down_10 = topGO(go, ont="BP", sort="Down", n=10)
    print(top_up_10)
    print(top_down_10)
    return(go)
}