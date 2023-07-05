#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# 06FunctionalAnalysis.R
#
#
# Functional Analysis:
#
#     This script allows to carry out a functional analysis based on 
#     two different methods: ORA and GSEA.
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

parser <- ArgumentParser(description="This script compute functional 
                                      analysis on meta-analysis results")

# Data directory
parser$add_argument("-i", "--indir",
                    action="store",
                    type="character",
                    default=".",
                    help="Data directory, by default the current
                          directory will be taken.")

# Output directory
parser$add_argument("-o", "--outdir",
                    action="store",
                    type="character",
                    default=".",
                    help="Where you would like the output files to be placed,
                          by default the current directory will be taken.")

# Output name
parser$add_argument("-n", "--name",
                    action="store",
                    type="character",
                    default="FA",
                    help="Output name")

# Input files path
parser$add_argument("-f", "--files",
                    action="store",
                    type="character",
                    required=TRUE,
                    help="Input files path")

# Contrasts
parser$add_argument("-c", "--contrasts",
                    action="store",
                    type="character",
                    required=TRUE,
                    help="Contrasts to include.") 

# Method
parser$add_argument("-m", "--method",
                    action="store",
                    type="character",
                    choices = c("ORA","GSEA"),
                    default= "GSEA",
                    help="Functional anlysis method")  

# Database 
parser$add_argument("-d", "--database",
                    action="store",
                    type="character",
                    default= "Reactome",
                    help="Data base to use. One or more of:
                          BP, MF, CC, KEGG, Reactome")  

# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Execution 
args <- parser$parse_args(args = c('-i=/home/rmoldovan/T2D-Meta-Analysis/Data/MA/Obesity',
                                   '-o=/home/rmoldovan/T2D-Meta-Analysis/Data/FA',
                                   '-f=c("IO/Meta-analysis_1_DF.RData",   \\
                                         "IOM/Meta-analysis_2_DF.RData",  \\
                                         "IOF/Meta-analysis_3_DF.RData",   \\
                                         "SDIO/Meta-analysis_4_DF.RData")',
                                   '-c=c("Ob_IS - Np_IS",  \\
                                         "Ob_IS.M - Np_IS.M",  \\
                                         "Ob_IS.F - Np_IS.F",  \\
                                         "(Ob_IS.M - Np_IS.M) - (Ob_IS.F-  Np_IS.F)")',
                                   '-d=c("BP", "MF", "CC", "KEGG", "Reactome")')) 

#------------- Checking arguments
args$files = eval(parse(text = args$files))
args$contrasts = eval(parse(text = args$contrasts))
args$database = eval(parse(text = args$database)) 

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
    Data[[args$isec[[1]]]]$Up = setdiff(Data[[args$isec[[1]]]][[1]],
                                        Data[[args$isec[[2]]]][[1]])
    Data[[args$isec[[1]]]]$Down = setdiff(Data[[args$isec[[1]]]][[2]],
                                          Data[[args$isec[[2]]]][[2]])
    Data[[args$isec[[2]]]]$Up = setdiff(Data[[args$isec[[2]]]][[1]],
                                        Data[[args$isec[[1]]]][[1]])
    Data[[args$isec[[2]]]]$Down = setdiff(Data[[args$isec[[2]]]][[2]],
                                          Data[[args$isec[[2]]]][[1]])
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
    
    save(Results, file=glue("{DirOut}/{args$name}_oraBP.RData"))
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
    
    save(Results, file=glue("{DirOut}/{args$name}_oraMF.RData"))
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
    
    save(Results, file=glue("{DirOut}/{args$name}_oraCC.RData"))
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
    
    save(Results, file=glue("{DirOut}/{args$name}_oraKEGG.RData"))
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
    
    save(Results, file=glue("{DirOut}/{args$name}_oraReactome.RData"))
  }
}else{
  cat("\n> GSEA\n")
  cat("\n---------------------------------------------\n")
  # Prepare the data  
  # Gene lists
  Data = sapply(Datas, function(x){
    geneList <- x$logFC
    #geneList <- sapply(seq(nrow(x)), function(i) ifelse( x[i, "logFC"] > 0, x[i, "p.adjust.BH"], -x[i, "p.adjust.BH"]))
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
    
    save(Results, file=glue("{DirOut}/{args$name}_gseBP.RData"))
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
    
    save(Results, file=glue("{DirOut}/{args$name}_gseMF.RData"))
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
    
    save(Results, file=glue("{DirOut}/{args$name}_gseCC.RData"))
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
    
    save(Results, file=glue("{DirOut}/{args$name}_gseKEGG.RData"))
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
    
    save(Results, file=glue("{DirOut}/{args$name}_gseReactome.RData"))
  }
  
}

cat("\nDone\n")
cat("\n---------------------------------------------\n")

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
