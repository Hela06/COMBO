#!/usr/bin/env Rscript

print("Load packages!")
suppressPackageStartupMessages({
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(plyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(biomaRt)
})
#make.input.neo4j <- function (edge.CH3, edge.mRNA, edge.meta, node.meta, gene.exp, CH3.status, lfm, lfe, adjpval, out.edge, out.node) {
  myargs = commandArgs(trailingOnly=TRUE)
  
  print("PRELIMINARY PHASE: ANNOTATION FILE CREATION")
  print("Load annotation Illumina 450K")
  ann450k_hg19 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k_hg19 <- data.frame(ann450k_hg19)
  ann450k_sub <- ann450k_hg19[,c(1,4,24,26)]
  ann450k_sub$chr <- gsub("chr", "", ann450k_sub$chr)
  colnames(ann450k_sub)[1] <- "CHR"
  
  
  print("Load CH3 status")
  condition <- myargs[1]
  path <- myargs[2]
  CH3.status <- paste0(path, "/DMProbes_", condition, "_limma.txt")
  CH3_status <- read.delim(CH3.status) 
  adjpval <- myargs[3]
  CH3_status <- subset(CH3_status, adj.P.Val < adjpval)
  CH3_status <- CH3_status[,c(4,9)]
  
  lfm <- myargs[4]
  lfm <- as.numeric(lfm)
  if(lfm == 0){
  CH3_status$CH3_STATE <- ifelse(CH3_status$logFC > 0, "HYPER", "HYPO")
  } else {
    CH3_status <- CH3_status(subset, logFC < - lfm | logFC > lfm )
    CH3_status$CH3_STATE <- ifelse(CH3_status$logFC > 0, "HYPER", "HYPO")
  }
  
  print("Create methylation annotation")
  ann450k_sub <- left_join(ann450k_sub, CH3_status, by="Name")
  ann450k_sub$logFC <- NULL
  ann450k_sub[is.na(ann450k_sub)] <- ""
  ann450k_sub$UCSC_RefGene_Name <- ifelse(ann450k_sub$UCSC_RefGene_Name == "", ann450k_sub$Name, ann450k_sub$UCSC_RefGene_Name)
  
  ann450k_sub <- separate_rows(ann450k_sub,c(3,4),sep = ";")
  colnames(ann450k_sub)[3:4] <- c("SYMBOL", "CH3_POSITION")
  
  ann450k_sub$HYPO <- ifelse(ann450k_sub$CH3_STATE == "HYPO", ann450k_sub$CH3_POSITION, "")
  ann450k_sub$HYPER <- ifelse(ann450k_sub$CH3_STATE == "HYPER", ann450k_sub$CH3_POSITION, "")
  
  
  ann450k_sub[,c(4,5)] <- NULL
  
  ann450k_sub$ENTREZID <- mapIds(org.Hs.eg.db, keys=ann450k_sub$SYMBOL, column="ENTREZID", keytype="SYMBOL", multiVals="first")
  ann450k_sub$ENTREZID <- ifelse(is.na(ann450k_sub$ENTREZID),ann450k_sub$SYMBOL, ann450k_sub$ENTREZID)
  ann450k_sub$SYMBOL <- as.character(ann450k_sub$SYMBOL)
  ann450k_sub <- ann450k_sub[order(ann450k_sub$Name),]
  ann450k_sub <- distinct(ann450k_sub)
  
  ann450k_two <- ddply(ann450k_sub, .(ENTREZID, SYMBOL, CHR), summarize,
                       HYPO=paste(HYPO[HYPO != ""],collapse="/"),
                       HYPER=paste(HYPER[HYPER != ""],collapse="/"))
                       
  probe <- as.vector(ann450k_sub$Name)
  ann450k_two <- ann450k_two[!(ann450k_two$SYMBOL %in% probe),]
  
  ann450k_sub <- ddply(ann450k_sub, .(ENTREZID, SYMBOL, CHR), summarize,
                       HYPO=paste(HYPO[HYPO != ""],collapse="/"),
                       HYPER=paste(HYPER[HYPER != ""],collapse="/"),
                       Name=paste(Name, collapse="/"))
  ann450k_sub$'SEQ:ID' <- seq(1:nrow(ann450k_sub))
  ann450k_sub <- separate_rows(ann450k_sub,6,sep = "/")
  #ann450k_sub <- distinct(ann450k_sub)
  
  #GENE EXPRESSION 
  print("Load annotation RNA_seq data")
  dataset="hsapiens_gene_ensembl"
  mart = biomaRt::useMart(biomart="ensembl", dataset = dataset, host="https://apr2022.archive.ensembl.org") 
  gene.annotations <- biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "chromosome_name"))
  colnames(gene.annotations)[1:4] <- c( "ENSEMBL","SYMBOL", "ENTREZID", "CHR")
  
  print("Load Gene expression")
  #controllare dettagliatamente le colonne
  gene.exp <- paste0(path, "/DEgenes_", condition, "_limma.txt")
  lfe <- myargs[5]
  lfe <- as.numeric(lfe)
  GENE_exp <- read.delim(gene.exp)
  GENE_exp <- GENE_exp[,c(1,2,6,8,9,10)]
  colnames(GENE_exp)[1] <- "ENSEMBL"
  GENE_exp$GENEXP <- ifelse(GENE_exp$logFC >= lfe & GENE_exp$adj.P.Val < adjpval, "UP_REG", "DOWN_REG")
  GENE_exp$GENEXP <- ifelse(between(GENE_exp$logFC, -lfe, lfe), "", GENE_exp$GENEXP)
  
  GENE_exp$adj.P.Val <- NULL
  GENE_exp$ENTREZID <- as.character(GENE_exp$ENTREZID)
  
  
  print("Create total annotation")
  TOT_ann <- left_join(ann450k_sub, GENE_exp, by=c("ENTREZID","SYMBOL", "CHR")) 
  TOT_ann <- distinct(TOT_ann)
  TOT_ann <- ddply(TOT_ann, .(ENTREZID, SYMBOL, CHR, HYPO, HYPER, Name, `SEQ:ID`, ENSEMBL), summarize,
                   logFC=paste(logFC,collapse="/"),
                   GENEXP=paste(GENEXP,collapse="/"))
                   
  TOT_ann <- distinct(TOT_ann)
  
  print("Create annotation for gene expression nodes")
  gene.annotations$'SEQ:ID' <- seq(from = (nrow(ann450k_sub)+1), length.out=nrow(gene.annotations))
  gene.annotations$ENTREZID <- as.character(gene.annotations$ENTREZID)
  TOT_ann_exp <- left_join( GENE_exp,ann450k_two, by=c("ENTREZID","SYMBOL", "CHR"))
  TOT_ann_exp <- distinct(TOT_ann_exp)
  
  TOT_ann_exp <- left_join( gene.annotations, TOT_ann_exp, by=c("ENSEMBL", "ENTREZID", "SYMBOL", "CHR"))
  TOT_ann_exp <- distinct(TOT_ann_exp)
  
  
  #loading custom annotation
  myarray <- unlist(strsplit(myargs[6],","))
  name_ann <- unlist(strsplit(myargs[7],","))
  if(myarray == "NULL"){
    print("No custom annotation to add!")
  } else {
  print("Load and add custom annotation")
  for (n in 1:length(myarray)) {
    df <- read.delim(myarray[n], na.strings = "")
    if(ncol(df) <= 2){
      colnames(df)[2] <- name_ann[n]
      TOT_ann <- left_join(TOT_ann,df, by=c("ENSEMBL")) 
      TOT_ann_exp <- left_join(TOT_ann_exp,df, by=c("ENSEMBL"))
    }
    else if (ncol(df) > 2){
      df <- df[,c("ENSEMBL", condition)]
      colnames(df)[2] <- name_ann[n]
      TOT_ann <- left_join(TOT_ann,df, by=c("ENSEMBL")) 
      TOT_ann_exp <- left_join(TOT_ann_exp,df, by=c("ENSEMBL"))
    }
  }
  }
  
  
  print("Refine final annotation")
  TOT_ann[is.na(TOT_ann)] <- ""
  TOT_ann[TOT_ann == "NA"] <- ""
  TOT_ann_exp[is.na(TOT_ann_exp)] <- ""
  TOT_ann_exp[TOT_ann_exp == "NA"] <- ""
  
  print("PHASE1: EDGEs CREATION")
  print("Process CH3_edge")
  path_b <- myargs[8]
  edge.CH3 <- paste0(path_b, "/Mval_",condition, "_BNout.txt")
  CH3_edge <- read.delim(edge.CH3)
  CH3_edge[,c(2,3)] <- NULL
  CH3_edge <- separate(CH3_edge, Implication, into = c(LETTERS[1:11]), sep = " ")
  CH3_edge[is.na(CH3_edge)] <- 0
  
  CH3_edge$F <- apply(CH3_edge, 1, function(row){
    if(row["F"] == "AND"){
      if(row["B"] == "low" & row["E"] == "low" & row["H"] == "high" & row["K"] == "high"){
        row["F"] <- "both_LL-HH"
      }
      else if(row["B"] == "low" & row["E"] == "high" & row["H"] == "high" & row["K"] == "low"){
        row["F"] <- "both_LH-HL"
      }
    }
    else{row["F"] <- paste(row["B"], row["E"], sep = "-")
    }
  })
  
  CH3_edge[,c(2,5,7:11)] <- NULL 
  CH3_edge[,c(2)] <- NULL
  CH3_edge <- CH3_edge[CH3_edge$D != 0, ]
  
  print("Process mRNA edge")
  edge.mRNA <- paste0(path_b, "/log2fpkm_",condition, "_BNout.txt")
  GE_edge <- read.delim(edge.mRNA)
  GE_edge[,c(2,3)] <- NULL
  GE_edge <- separate(GE_edge, Implication, into = c(LETTERS[1:11]), sep = " ")
  GE_edge[is.na(GE_edge)] <- 0
  GE_edge$F <- apply(GE_edge, 1, function(row){
    if(row["F"] == "AND"){
      if(row["B"] == "low" & row["E"] == "low" & row["H"] == "high" & row["K"] == "high"){
        row["F"] <- "both_LL-HH"
      }
      else if(row["B"] == "low" & row["E"] == "high" & row["H"] == "high" & row["K"] == "low"){
        row["F"] <- "both_LH-HL"
      }
    }
    else{row["F"] <- paste(row["B"], row["E"], sep = "-")
    }
  })
  
  GE_edge[,c(2,5,7:11)] <- NULL 
  GE_edge[,c(2)] <- NULL
  GE_edge <- GE_edge[GE_edge$D != 0, ]
  
  print("PHASE2: NODEs CREATION")
  print("Process node CH3 and mRNA")
  node_CH3 <- c(CH3_edge$A, CH3_edge$D)
  node_CH3 <- node_CH3[!duplicated(node_CH3)]
  node_GE <- c(GE_edge$A, GE_edge$D)
  node_GE <- node_GE[!duplicated(node_GE)]
  node_CH3 <- data.frame(node_CH3)
  node_CH3$"label:LABEL" <- "CH3"
  
  node_CH3 <- left_join(node_CH3, TOT_ann, by=c("node_CH3" = "Name"))
  #node_CH3 <- node_CH3[!duplicated(node_CH3),]
  
  #node_CH3$node_CH3 <- NULL
  
  node_GE <- data.frame(node_GE)
  node_GE$"label:LABEL" <- "mRNA"
  colnames(node_GE)[1] <- "ENSEMBL"
  node_GE$ENSEMBL <- as.character(node_GE$ENSEMBL)
  node_GE <- left_join(node_GE, TOT_ann_exp, by="ENSEMBL")
  
  
  print("Load meta node and merge")
  meta <- myargs[9]
  path_m <- myargs[10]
  if(meta == "KEGG"){
    node.meta <- paste0(path_m, "/KEGG/nodes.txt")
    edge.meta <- paste0(path_m, "/KEGG/edges.txt")
    node_meta <- read.delim(node.meta)
    edge_meta <- read.delim(edge.meta)
  } else if (meta == "KEGG_Reactome"){
    node.meta <- paste0(path_m, "/KEGG_Reactome/nodes.txt")
    edge.meta <- paste0(path_m, "/KEGG_Reactome/edges.txt")
    node_meta <- read.delim(node.meta)
    edge_meta <- read.delim(edge.meta)
  }
    
  node_meta <- as.data.frame(node_meta)
  node_meta$"label:LABEL" <- "MP"
  node_meta$X.Id <- as.character(node_meta$X.Id)
  node_meta <- left_join(node_meta, TOT_ann_exp, by=c("X.Id" = "ENTREZID"))
  node_meta$'SEQ:ID' <- seq(from = (nrow(ann450k_sub)+ nrow(gene.annotations)+1), length.out=nrow(node_meta))
  
  node_meta$`SEQ:ID` <- ifelse(node_meta$Type != "GENE", node_meta$X.Id, node_meta$`SEQ:ID`) 
  node_meta[,c(3:4)] <- NULL
  colnames(node_meta)[1:2] <- c("ENTREZID", "SYMBOL")
  node_meta[,5] <- NULL
  
  print("Insert unique ID on edges")
  
  IDENTIFIERS_CH3 <- node_CH3[,c("node_CH3","SEQ:ID")]
  node_CH3$node_CH3 <- NULL
  CH3_edge <- left_join(CH3_edge, IDENTIFIERS_CH3, by=c("A" = "node_CH3"))
  CH3_edge <- left_join(CH3_edge, IDENTIFIERS_CH3, by=c("D" = "node_CH3"))
  CH3_edge[,c(1,2)] <- NULL
  colnames(CH3_edge) <- c("ATTR:TYPE", "SRG:START_ID", "TRG:END_ID")
  CH3_edge$ATTR2 <- NA
  
  IDENTIFIERS_GE <- node_GE[,c("ENSEMBL","SEQ:ID")]
  GE_edge <- left_join(GE_edge, IDENTIFIERS_GE, by=c("A" = "ENSEMBL"))
  GE_edge <- left_join(GE_edge, IDENTIFIERS_GE, by=c("D" = "ENSEMBL"))
  GE_edge[,c(1,2)] <- NULL
  colnames(GE_edge) <- c("ATTR:TYPE", "SRG:START_ID", "TRG:END_ID")
  GE_edge$ATTR2 <- NA
  node_GE$ENTREZID <- ifelse(is.na(node_GE$ENTREZID) | node_GE$ENTREZID == "", node_GE$SYMBOL, node_GE$ENTREZID)
  
  
  IDENTIFIERS_META <- node_meta[,c("ENTREZID","SEQ:ID")]
  edge_meta[,5] <- NULL
  edge_meta <- left_join(edge_meta, IDENTIFIERS_META, by= c("X.Start" = "ENTREZID"))
  edge_meta <- left_join(edge_meta, IDENTIFIERS_META, by= c( "End" = "ENTREZID"))
  edge_meta[,c(1,2)] <- NULL
  colnames(edge_meta) <- c( "ATTR2","ATTR:TYPE","SRG:START_ID", "TRG:END_ID")
  #node_meta <- node_meta[ ! node_meta$ENTREZID %in% entrez_to_remove,] 
  
  print("Insert interlayer edges")
  
  IDENTIFIERS_COMMON <- node_CH3[,c("ENSEMBL","SEQ:ID")]
  INTERL_GC <- inner_join(IDENTIFIERS_GE, IDENTIFIERS_COMMON, by="ENSEMBL")
  colnames(INTERL_GC) <- c("ATTR:TYPE", "SRG:START_ID", "TRG:END_ID")
  INTERL_GC$`ATTR:TYPE` <- "INTERLAYER"
  INTERL_CG <- inner_join(IDENTIFIERS_COMMON, IDENTIFIERS_GE, by="ENSEMBL")
  colnames(INTERL_CG) <- c("ATTR:TYPE", "SRG:START_ID", "TRG:END_ID")
  INTERL_CG$`ATTR:TYPE` <- "INTERLAYER"
  
  IDENTIFIERS_COMMON <- node_CH3[,c("ENTREZID","SEQ:ID")]
  IDENTIFIERS_GE <- node_GE[,c("ENTREZID","SEQ:ID")]
  INTERL_GM <- inner_join(IDENTIFIERS_GE, IDENTIFIERS_META, by="ENTREZID")
  colnames(INTERL_GM) <- c("ATTR:TYPE", "SRG:START_ID", "TRG:END_ID")
  INTERL_GM$`ATTR:TYPE` <- "INTERLAYER"
  INTERL_MG <- inner_join(IDENTIFIERS_META, IDENTIFIERS_GE, by="ENTREZID")
  colnames(INTERL_MG) <- c("ATTR:TYPE", "SRG:START_ID", "TRG:END_ID")
  INTERL_MG$`ATTR:TYPE` <- "INTERLAYER"
  
  INTERL_CM <- inner_join(IDENTIFIERS_COMMON, IDENTIFIERS_META, by="ENTREZID")
  colnames(INTERL_CM) <- c("ATTR:TYPE", "SRG:START_ID", "TRG:END_ID")
  INTERL_CM$`ATTR:TYPE` <- "INTERLAYER"
  INTERL_MC <- inner_join(IDENTIFIERS_META, IDENTIFIERS_COMMON, by="ENTREZID")
  colnames(INTERL_MC) <- c("ATTR:TYPE", "SRG:START_ID", "TRG:END_ID")
  INTERL_MC$`ATTR:TYPE` <- "INTERLAYER"
  
  INTERLAYER_TOT <- rbind(INTERL_CG, INTERL_GC, INTERL_CM, INTERL_MC, INTERL_GM, INTERL_MG)
  INTERLAYER_TOT$ATTR2 <- NA
  INTERLAYER_TOT <- distinct(INTERLAYER_TOT)
  
  node_GE$ENSEMBL <- NULL
  node_CH3$ENSEMBL <- NULL
  node_meta$ENSEMBL <- NULL
  
  print("Merge edge")
  edge_tot <- rbind(CH3_edge, GE_edge, edge_meta, INTERLAYER_TOT)
  edge_tot <- distinct(edge_tot)
  edge_tot <- edge_tot[edge_tot$`SRG:START_ID` != edge_tot$`TRG:END_ID`,]
  print("Save edge ")
  path_out <- myargs[11]
  out.edge <- paste0(path_out, "/", condition, "_multilayer_edges.csv")
  write.csv(edge_tot, file = out.edge, quote = FALSE, sep = ",", row.names = FALSE)
  
  print("Merge node")
  #write.csv(node_CH3, file = "./node_CH3.csv", quote = FALSE, sep = ",", row.names = FALSE)
  #write.csv(node_GE, file = "./node_GE.csv", quote = FALSE, sep = ",", row.names = FALSE)
  #write.csv(node_meta, file = "./node_meta.csv", quote = FALSE, sep = ",", row.names = FALSE)
  
  node_tot <- rbind(node_GE, node_CH3 ,node_meta)
  node_tot <- distinct(node_tot)
  
  #node_tot <- ddply(node_tot, .(`SEQ:ID`,ENTREZID, SYMBOL,`label:LABEL`, MUTATION, PROTEXP, CH3_POSITION, CH3_STATE, CH3_STATE_GvD, logFC, logFC_GvD, GENEXP, GENEXP_GvD, TRANFAC), summarize,
  #MAP=paste(MAP,collapse="/"))
  
  #rimuovere duplicati lasciando logFC pi? alti
  print("Remove duplicated")
  node_tot$logFC <- as.numeric(node_tot$logFC)
  node_tot <- node_tot[order(node_tot$`SEQ:ID`, -abs(node_tot$logFC) ), ]
  node_tot <- node_tot[ !duplicated(node_tot$`SEQ:ID`), ]
  
  node_tot$SYMBOL <- gsub(",", "/", node_tot$SYMBOL)
  node_tot[is.na(node_tot)] <- "" 
  node_tot[node_tot == "NA"] <- "" 
  #node_tot$SYMBOL <- as.character(node_tot$SYMBOL)
  #node_tot$ENTREZID <- ifelse(is.na(node_tot$ENTREZID), node_tot$SYMBOL, node_tot$ENTREZID)
  
  print("Save node")
  out.node <- paste0(path_out, "/", condition, "_multilayer_nodes.csv")
  #write.csv(node_neo, file = "./Neo4J/node_tot_neo4j.csv", quote = FALSE, sep = ",", row.names = FALSE)
  write.csv(node_tot, file = out.node, quote = FALSE, sep = ",", row.names = FALSE)
#}


