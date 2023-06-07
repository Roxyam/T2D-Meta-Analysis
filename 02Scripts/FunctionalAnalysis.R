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
set.seed(1808)

# ~~~~~~~~~~~~ Parameters ~~~~~~~~~~~~ #

# <- ArgumentParser(description="")

# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Prepare arguments
args = list()
args$indir = "C:/Users/roxya/OneDrive/Documentos/01Master_bioinformatica/00TFM/Git/T2D-Meta-Analysis/Data/MA/Obesity"
#args$indir = "C:/Users/roxya/OneDrive/Documentos/01Master_bioinformatica/00TFM/Git/T2D-Meta-Analysis/Data/MA/Diabetes"

args$outdir = "C:/Users/roxya/OneDrive/Documentos/01Master_bioinformatica/00TFM/Git/T2D-Meta-Analysis/Data/FA_intersect"

args$outname = "Obesidad"
#args$outname = "Diabetes"

args$files = c("ObvsNp/Meta-analysis_1_DF.RData", 
                  "ObMvsNpM/Meta-analysis_2_DF.RData",
                  "ObFvsNpF/Meta-analysis_3_DF.RData",
                  "ObM_NpMvsObF_NpF/Meta-analysis_4_DF.RData")
'args$files = c("ObM_NpMvsObF_NpF/Meta-analysis_4_DF.RData")
args$files = c("ObIRvsObIS/Meta-analysis_1_DF.RData", 
               "ObIRMvsObISM/Meta-analysis_2_DF.RData",
               "ObIRFvsObISF/Meta-analysis_3_DF.RData",
               "IRM_ISMvsIRF_ISF/Meta-analysis_4_DF.RData")'

args$contrasts = c("Ob_IS - Np_IS", "Ob_IS.M - Np_IS.M", "Ob_IS.F - Np_IS.F", "(Ob_IS.M - Np_IS.M) - (Ob_IS.F-  Np_IS.F)")
#args$contrasts = c( "(Ob_IS.M - Np_IS.M) - (Ob_IS.F-  Np_IS.F)")
'args$contrasts = c("Ob_IR - Ob_IS", "Ob_IR.M - Ob_IS.M", "Ob_IR.F - Ob_IS.F",
                   "(Ob_IR.M-Ob_IS.M) - (Ob_IR.F-Ob_IS.F)")'

args$method = c("GSEA") # ORA or GSEA
args$database = c("BP", "MF", "CC", "KEGG", "Reactome")
args$isec = c(2,3) # Si ponemos intersect en los de sexo toma solo los exclusivos?

# Get data
DirOut = args$outdir
files = glue("{args$indir}/{args$files}")
Datas = sapply(files, function(x) get(load(file=x)), simplify = FALSE)
contrasts = args$contrasts

#------------- Functional Analysis
if(args$method == "ORA"){
  cat("\n> ORA\n")
  cat("\n---------------------------------------------\n")
  # Prepare the data
  p_lim = 0.05
  logFC_lim = 0
  
  # Gene lists 
  
  Data = sapply(Datas, function(Data){
    Up = rownames(Data[which(Data$p.adjust.BH < p_lim & Data$logFC > logFC_lim), ])
    Down = rownames(Data[which(Data$p.adjust.BH < p_lim & Data$logFC < 0), ])
    #Up = (Data[which(Data$p.adjust.BH < p_lim & Data$logFC > logFC_lim), "Symbol"])
    #Down = (Data[which(Data$p.adjust.BH < p_lim & Data$logFC < 0), "Symbol"])
    return(c("Up" = list(Up), "Down" = list(Down)))
  }, simplify = FALSE)
  names(Data) = contrasts
  
  if(! is.null(args$isec)){
    Data[[args$isec[[1]]]]$Up = setdiff(Data[[args$isec[[1]]]][[1]], Data[[args$isec[[2]]]][[1]])
    Data[[args$isec[[1]]]]$Down = setdiff(Data[[args$isec[[1]]]][[2]], Data[[args$isec[[2]]]][[2]])
    Data[[args$isec[[2]]]]$Up = setdiff(Data[[args$isec[[2]]]][[1]], Data[[args$isec[[1]]]][[1]])
    Data[[args$isec[[2]]]]$Down = setdiff(Data[[args$isec[[2]]]][[2]], Data[[args$isec[[2]]]][[1]])
      }
  
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
  cat("\n> GSEA\n")
  cat("\n---------------------------------------------\n")
  # Prepare the data  
  # Gene lists
  Data = sapply(Datas, function(x){
    geneList <- x$logFC
    names(geneList) <- as.character(row.names(x))
    geneList <- sort(geneList, decreasing = TRUE) 
    return(geneList)
  }, simplify = FALSE)   
  names(Data) = contrasts
  
  # Go BP
  if ("BP" %in% args$database){ 
    cat("\n\tBP\n") 
    Results = mclapply(contrasts, function(c){
      geneList = Data[[c]]
      enrich_result <- gseGO(geneList      = geneList,
                             OrgDb         = org.Hs.eg.db,  
                             keyType       = "ENTREZID",  
                             ont           = "BP", 
                             pAdjustMethod = "BH",  
                             pvalueCutoff  = 0.05,
                             seed = TRUE)
      return(enrich_result)})
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_gseBP.RData"))
  }
  # GO MF
  if ("MF" %in% args$database){ 
    cat("\n\tMF\n") 
    Results = mclapply(contrasts, function(c){
      geneList = Data[[c]]
      enrich_result <- gseGO(geneList      = geneList,
                             OrgDb         = org.Hs.eg.db,  
                             keyType       = "ENTREZID",  
                             ont           = "MF", 
                             pAdjustMethod = "BH",  
                             pvalueCutoff  = 0.05,
                             seed = TRUE)
      return(enrich_result)})
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_gseMF.RData"))
  }
  # GO CC
  if ("CC" %in% args$database){ 
    cat("\n\tCC\n") 
    Results = mclapply(contrasts, function(c){
      geneList = Data[[c]]
      enrich_result <- gseGO(geneList      = geneList,
                             OrgDb         = org.Hs.eg.db,  
                             keyType       = "ENTREZID",  
                             ont           = "CC", 
                             pAdjustMethod = "BH",  
                             pvalueCutoff  = 0.05,
                             seed = TRUE)
      return(enrich_result)}) 
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_gseCC.RData"))
  }
  # KEGG
  if ("KEGG" %in% args$database){ 
    cat("\n\tKEGG\n") 
    Results = mclapply(contrasts, function(c){
      geneList = Data[[c]]
      enrich_result <- gseKEGG(geneList      = geneList,
                               organism      = "hsa",
                               keyType       = "ncbi-geneid",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               seed = TRUE,
                               use_internal_data = TRUE)
      return(enrich_result)}) 
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_gseKEGG.RData"))
  }
  # Reactome
  if ("Reactome" %in% args$database){ 
    cat("\n\tReactome\n") 
    Results = mclapply(contrasts, function(c){
      geneList = Data[[c]]
      enrich_result <- gsePathway(geneList      = geneList,
                                  organism      = "human", 
                                  pAdjustMethod = "BH",  
                                  pvalueCutoff  = 0.05,
                                  seed = TRUE)
      return(enrich_result)})
    Results = setNames(Results, contrasts) 
    
    save(Results, file=glue("{DirOut}/{args$outname}_gseReactome.RData"))
  }
  
}



