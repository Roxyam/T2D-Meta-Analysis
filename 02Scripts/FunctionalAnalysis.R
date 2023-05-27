#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# 07FunctionalAnalysis.R
#
#
# Compute ORA or GSEA on Mata-Analysis results:
# 
#
#
#     ****************************************************************
#     *    Author:  Roxana Andreea Moldovan                          *
#     *    Contact: roxana.andreea.moldovan@gmail.com                *
#     *    Version: 1.0                                              *
#     *    Creation date: 27/05/2023                                 *
#     ****************************************************************
#
#
# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #


# ~~~~~~~~~~~~ Loading packages ~~~~~~~~~~~~ #

library(pacman)
pacman::p_load(argparse)
pacman::p_load(clusterProfiler)
pacman::p_load(dplyr)
pacman::p_load(glue)
pacman::p_load(org.Hs.eg.db)
pacman::p_load(parallel)
pacman::p_load(ReactomePA)
pacman::p_load(stringr)
pacman::p_load(sf)


# ~~~~~~~~~~~~ Functions ~~~~~~~~~~~~ #


# ~~~~~~~~~~~~ Parameters ~~~~~~~~~~~~ #

parser <- ArgumentParser(description="")

# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Prepare arguments
args = list()
args$indir = "C:/Users/roxya/OneDrive/Documentos/01Master_bioinformatica/00TFM/Git/T2D-Meta-Analysis/Data/MA/Obesity"
args$outdir = "C:/Users/roxya/OneDrive/Documentos/01Master_bioinformatica/00TFM/Git/T2D-Meta-Analysis/Data/FA"
args$outname = "Obesidad"
args$files = c("ObvsNp/Meta-analysis_1_DF.RData", 
                  "ObMvsNpM/Meta-analysis_2_DF.RData",
                  "ObFvsNpF/Meta-analysis_3_DF.RData",
                  "ObM_NpMvsObF_NpF/Meta-analysis_4_DF.RData")
args$contrasts = c("Ob - Np", "Ob.M - Np.M", "Ob.F - Np.F",
                   "(Ob.M - Np.M) - (Ob.F-  Np.F)")
args$method = c("ORA") # ORA or GSEA
args$database = c("BP", "MF", "CC", "KEGG", "Reactome")

# Get data
DirOut = args$outdir
files = glue("{args$indir}/{args$files}")
Datas = sapply(files, function(x) get(load(file=x)), simplify = FALSE)
contrasts = args$contrasts

#------------- Functional Analysis
if(args$method == "ORA"){
  cat("\n\t> ORA\n")
  cat("\n---------------------------------------------\n")
  # Prepare the data
  p_lim = 0.05
  logFC_lim = 0
  
  # Gene lists
  Direction = c("Up", "Down")
  
  Data = sapply(Datas, function(Data){
    Up = rownames(Data[which(Data$p.adjust.BH < p_lim & Data$logFC > logFC_lim), ])
    Down = rownames(Data[which(Data$p.adjust.BH < p_lim & Data$logFC < 0), ])
    return(c("Up" = list(Up), "Down" = list(Down)))
  }, simplify = FALSE)
  names(Data) = contrasts
  
  # Go BP
  if ("BP" %in% args$database){ 
    cat("\n\tBP\n") 
    Results = mclapply(contrasts, function(c){
      r = sapply(Data[[c]], function(x){
        enrich_result <- enrichGO(gene          = x,
                                  OrgDb         = org.Hs.eg.db,  
                                  keyType       = "ENTREZID",  
                                  ont           = "BP", 
                                  pAdjustMethod = "BH",  
                                  pvalueCutoff  = 0.05,  
                                  qvalueCutoff  = 0.05)
        return(enrich_result)}, simplify = FALSE, USE.NAMES = TRUE)
      return(r)
    })
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_oraBP.RData"))
  }
  # GO MF
  if ("MF" %in% args$database){ 
    cat("\n\tMF\n") 
    Results = mclapply(contrasts, function(c){
      r = sapply(Data[[c]], function(x){
        enrich_result <- enrichGO(gene          = x,
                                  OrgDb         = org.Hs.eg.db,  
                                  keyType       = "ENTREZID",  
                                  ont           = "MF", 
                                  pAdjustMethod = "BH",  
                                  pvalueCutoff  = 0.05,  
                                  qvalueCutoff  = 0.05)
        return(enrich_result)}, simplify = FALSE, USE.NAMES = TRUE)
      return(r)
    })
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_oraMF.RData"))
  }
  # GO CC
  if ("CC" %in% args$database){ 
    cat("\n\tCC\n") 
    # Realiza el análisis de enriquecimiento funcional
    Results = mclapply(contrasts, function(c){
      r = sapply(Data[[c]], function(x){
        enrich_result <- enrichGO(gene          = x,
                                  OrgDb         = org.Hs.eg.db,  
                                  keyType       = "ENTREZID",  
                                  ont           = "CC", 
                                  pAdjustMethod = "BH",  
                                  pvalueCutoff  = 0.05,  
                                  qvalueCutoff  = 0.05)
        return(enrich_result)}, simplify = FALSE, USE.NAMES = TRUE)
      return(r)
    })
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_oraCC.RData"))
  }
  # KEGG
  if ("KEGG" %in% args$database){ 
    cat("\n\tKEGG\n") 
    Results = mclapply(contrasts, function(c){
      r = sapply(Data[[c]], function(x){
        enrich_result <- clusterProfiler::enrichKEGG(gene          = x,
                                                     organism      = "hsa",  
                                                     keyType       = "ncbi-geneid",
                                                     pAdjustMethod = "BH",
                                                     pvalueCutoff = 0.05,
                                                     qvalueCutoff = 0.5,
                                                     use_internal_data = TRUE)
        return(enrich_result)}, simplify = FALSE, USE.NAMES = TRUE)
      return(r)
    })
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_oraKEGG.RData"))
  }
  # Reactome
  if ("Reactome" %in% args$database){ 
    cat("\n\tReactome\n") 
    Results = mclapply(contrasts, function(c){
      r = sapply(Data[[c]], function(x){
        enrich_result <- enrichPathway(gene      = x,
                                       organism      = "human",  
                                       pAdjustMethod = "BH",  
                                       pvalueCutoff  = 0.05,  
                                       qvalueCutoff  = 0.05)
        return(enrich_result)}, simplify = FALSE, USE.NAMES = TRUE)
      return(r)
    })
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_oraReactome.RData"))
  }
}else{
  
}



