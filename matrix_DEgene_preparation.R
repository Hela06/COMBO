#!/usr/bin/env Rscript

myargs = commandArgs(trailingOnly=TRUE)
#input.preparation <- function(expression,methylation, out_directory, info, comparison){
  #input.preparation <- function(){
  print("Load packages!")
  suppressPackageStartupMessages({
  library(data.table)
  library(DGEobj.utils) 
  library("biomaRt")
  library("dplyr")
  library(lumi) #convert B in Mvalue
  library(limma)
  library(edgeR)
  library(tibble)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  })
  
  print("Load sample information!") 
  samples <- read.table(myargs[4], check.names = FALSE)
  groups <- read.table(myargs[5], check.names = FALSE)
  print("Load expression matrix!")
  exp_matrix <- read.table(myargs[1], check.names = FALSE)
  print("Load methylation matrix!")
  met_matrix <- read.table(myargs[2], check.names = FALSE)
  directory <- myargs[3]
  
  exp_matrix <- exp_matrix[, row.names(samples)]
  met_matrix <- met_matrix[, row.names(samples)]
  
  print("Identify differential expressed genes")
  d0 <- DGEList(exp_matrix)
  keep <- filterByExpr(d0, design = samples)
  d0 <- d0[keep,keep.lib.sizes=FALSE]
  d0 <- calcNormFactors(d0, method="TMM")
  cnt <- d0$counts
  y <- voom(d0, samples, plot = F)
  fit_exp <- lmFit(y, samples)
  #head(coef(fit_exp))
  contMatrix <- as.matrix(groups)
  data.fit.con <- contrasts.fit(fit_exp, contMatrix)
  data.fit.eb <- eBayes(data.fit.con)
  
  dataset="hsapiens_gene_ensembl"
  mart = biomaRt::useMart(biomart="ensembl", dataset = dataset, host="https://apr2022.archive.ensembl.org") 
  gene.annotations <- biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "chromosome_name","start_position", "end_position"))
  gene.annotations <- dplyr::transmute(gene.annotations, external_gene_name,  ensembl_gene_id, entrezgene_id, chromosome_name, length = end_position - start_position)
  colnames(gene.annotations)[1:4] <- c( "SYMBOL", "ENSEMBL", "ENTREZID", "CHR")
  gene.anotations.sb <- gene.annotations[,1:4]
  #gene.anotations.sb$ENSEMBL <- as.character(gene.anotations.sb$ENSEMBL)
  
  for (c in 1:ncol(groups)) {
    nm <- row.names(groups)[which(groups[,c] == 1)]
    EXP <- topTable(data.fit.eb, coef = c, number = Inf, adjust.method = "BH")
    EXP <- dplyr::left_join(rownames_to_column(EXP), gene.anotations.sb, by=c("rowname" = "ENSEMBL"))
    name <- paste0(directory,"/DEgenes_", nm,"_limma.txt")
    write.table(EXP, file = name, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
  }
  
  print("Identify CpG differentially methylated")
  #eliminare +Inf e -Inf
  met_matrix[met_matrix == 1] <- 0.999999
  met_matrix[met_matrix == 0] <- 0.000001
  met_matrix_M <- beta2m(met_matrix)
  fit <- lmFit(met_matrix_M, samples)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  
  ann450_hg19 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450_hg19 <- as.data.frame(ann450_hg19)
  ann450_hg19_sub <- ann450_hg19[,c(1:4,18,19,24,26)]
  
  for (c in 1:ncol(groups)) {
    nm <- row.names(groups)[which(groups[,c] == 1)]
    DMPs <- topTable(fit2, num=Inf, coef=c, genelist = ann450_hg19_sub)
    name <- paste0(directory,"/DMProbes_", nm,"_limma.txt")
    write.table(DMPs, file = name, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
  }
  
  
  print("Identify comparison groups...")
  for(c in colnames(groups)){
    print(paste0("Creation input file for ", c, " group"))
    case <- row.names(groups)[which(groups[,c] == 1)]
    ctrl <- row.names(groups)[which(groups[,c] == -1)]
    case_n <- c()
    ctrl_n <- c()
    
    for (s in colnames(samples)) { 
      for (r in 1:nrow(samples)) {
      
        if( s == ctrl & samples[r,s] == 1 ){
          ctrl_n <- c(ctrl_n, row.names(samples[r,]))
        }
        else if( s == case & samples[r,s] == 1 ){
          case_n <- c(case_n, row.names(samples[r,]))
        }
      }
    }
    #ctrl_exp <- subset(exp_matrix, colnames(exp_matrix) %in% ctrl_n)
    #Modifiche Gianni
    case_exp <- subset(cnt, colnames(cnt) %in% case_n)
    
    #ctrl_met <- subset(met_matrix, colnames(met_matrix) %in% ctrl_n)
    case_met <- subset(met_matrix_M, colnames(met_matrix) %in% case_n)
    
      gene.annotations_n <- gene.annotations %>% dplyr::filter(ENSEMBL %in% rownames(case_exp))
      diff <- setdiff(row.names(case_exp), gene.annotations_n$ENSEMBL)
      gene.annotations_n <- gene.annotations_n %>% distinct(ENSEMBL, .keep_all = TRUE)
      case_exp <- as.matrix(case_exp)
      case_exp <- case_exp[!(row.names(case_exp) %in% diff), ]
      gene.annotations_n <- gene.annotations_n[order(match(gene.annotations_n$ENSEMBL, rownames(case_exp))),]
      featureLength <- gene.annotations_n$length
      
      Log2FPKM <- convertCounts(case_exp,
                                unit = "fpkm",
                                geneLength = featureLength,
                                log = TRUE) 
      
      file_e <- paste0(directory, "/log2fpkm_", case,"_SMinput.txt")
      write.table(Log2FPKM, file = file_e, quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
      
      file_m <- paste0(directory,"/Mval_", case,"_SMinput.txt")
      write.table(case_met, file = file_m, quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
    
  
  }
  
#}

#samples <- read.table("info_lung.txt")
#groups <- read.table("comparison_lung.txt", check.names = FALSE)
#exp_matrix <- read.table("expression_matrix_example.txt")
#met_matrix <- read.table("methylation_matrix_sb_example.txt")
