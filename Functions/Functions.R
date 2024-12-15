#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# Functions.R
#
#
#     This file contains a set of functions that are relevant for the 
#     entire project development. These functions are organized based 
#     on the steps in which they have been utilized. Within each
#     function, you can find detailed information about its
#     functionality, parameters, and dependencies.
#
#
#     ****************************************************************
#     *    Author:  Roxana Andreea Moldovan                          *
#     *    Contact: roxana.andreea.moldovan@gmail.com                *
#     *    Version: 1.1                                              *
#     *    Creation date: 26/03/2023                                 *
#     ****************************************************************
#
#
# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #

# ~~~~~~~~~~~~ General ~~~~~~~~~~~~ #

#------------- RX_coind
RX_coind <- function(dict, 
                     lst, 
                     names = FALSE){
  # ****************************************************** 
  # Function: RX_coind()       
  #     
  #   This function takes a list and modifies it based on information
  #   from a dictionary.
  #
  #   Input:
  #   -----
  #    - dict (vector with names):
  #         Dictionary where the keys correspond to elements in 
  #         the `lst` list, and the values represent the desired
  #         modifications for those elements.
  #    - lst (list):
  #         List to be modified according to the dictionary's information.
  #    - names (logical, optional, default is False):  
  #         Determines whether the original names of the list elements 
  #         will be preserved.
  # 
  #   Output:
  #   -----
  #   - (vector) List modified according to the dictionary's information.
  # 
  # Dependencies:
  # -----
  # none
  #
  # ******************************************************
  
  lst = as.vector(lst)
  out = (sapply(lst, function(x) dict[[x]], USE.NAMES = names))
  return(out)
}

#------------- RX_ifelse
RX_ifelse <- function(cond, 
                      true, 
                      false){
  # ****************************************************** 
  # Function: RX_ifelse()       
  #     
  #   This function implements a simplified version of the ifelse
  #   function in R.
  #
  #   Input:
  #   -----
  #    - cond: 
  #         Logical condition that determines whether to return 
  #         the true or false result.
  #    - true: 
  #         Value to return if the condition is true.
  #    - false:
  #         Value to return if the condition is false.
  # 
  #   Output:
  #   -----
  #   - (logical) Value corresponding to the evaluated condition.
  # 
  # Dependencies:
  # -----
  # none
  #
  # ******************************************************
  
  
  if (cond){
    out = true
  }else{
    out = false
  }
  return(out)
}


# ~~~~~~~~~~~~ 01Download ~~~~~~~~~~~~ #

#------------- RX_GetDataGEO
RX_GetDataGEO <- function(StudyAcc,
                          dir=getwd()){
  # ****************************************************** 
  # Function: RX_GetDataGEO
  #     
  #   This function is used to retrieve data from the Gene Expression Omnibus
  #   (GEO) using the provided GEO accession. It downloads the necessary files,
  #   decompresses them, and saves the study information in a metadata file.
  #
  # Input:
  # -----
  #   - StudyAcc (character): 
  #         The GEO accession identifier of the study.
  #   - dir (character): 
  #         Working directory. By default, the current directory is used.
  # Output:
  # -----
  #   none
  #
  # Dependencies:
  # -----
  #   library(GEOquery)
  #
  # ******************************************************
  
  setwd(dir) 
  print(glue("Processing {StudyAcc}:"))
  
  # Check if the file exists
  if(
    class(
      try(
        {
          getGEOSuppFiles(StudyAcc, makeDirectory = TRUE,
                          baseDir = dir, fetch_files = TRUE)
          # Set the working directory
          setwd(glue("{dir}/{StudyAcc}"))
        }
      )) == "try-error"
  ){
    # Error
    print(glue("{StudyAcc} could not be downloaded."))
  }else{
    # Descomprimir archivos
    files = list.files()
    for (file in files){
      if (endsWith(file, "tar")){
        system(glue("tar -xvf {file}"))
      } else{
        if (endsWith(file, "tar.gz")){
          system(glue("tar -xzvf {file}"))
        } else {
          if (endsWith(file, "gz")){ 
            system(glue("gzip -d {file}"))
          }
        }
      }
    }
    # Delete compressed files
    system(paste("rm", files, sep = ' '))
    
    # Download study information
    dir.create("01RawData")
    metadata = getGEO(GEO = StudyAcc)[[1]]
    save(metadata, file = glue("{StudyAcc}_metadata.RData"))
    system( "mv `ls | grep -v 01RawData` 01RawData/")
    
    # Completed
    print(glue("Study {StudyAcc} successfully processed."))
  }
}

  
#------------- RX_GetDataArrEx
RX_GetDataArrEx <- function(StudyAcc, dir = getwd()){
  # ****************************************************** 
  # Function: RX_GetDataArrEx
  #     
  #   This function is used to retrieve data from ArrayExpress using the
  #   provided ArrayExpress accession. It creates a directory to store the
  #   data, downloads the dataset, and saves the ExpressionSet object in an
  #   .rda file.
  #
  # Input:
  # -----
  #   - StudyAcc (character): 
  #         The ArrayExpress accession identifier of the study.
  #   - dir (character): 
  #         Working directory. By default, the current directory is used.
  # Output:
  # -----
  #   none
  #
  # Dependencies:
  # -----
  #   library(ArrayExpress)
  #
  # ******************************************************
  
  setwd(dir)
  # Create the directory to store the data
  path0 = glue("{StudyAcc}/01RawData")
  system(glue("mkdir {StudyAcc} {path0}"))
  
  # Data download
  if(
    class(
      try(
        {exp_set = ArrayExpress(StudyAcc, path = path0)} 
      )) == "try-error"){
    # Error
    print(glue("{StudyAcc} could not be downloaded."))
  }else{
    # Guardar el ExpressionSet
    save(exp_set, file = glue("{path0}/{StudyAcc}_raw.rda"))
  }
}

# ~~~~~~~~~~~~ 02Preprocessing ~~~~~~~~~~~~ #

#------------- RX_probe2gene
RX_probe2gene <- function(Data, 
                          anot_tb = NULL,
                          anotdb=org.Hs.eg.db,  
                          from = "PROBEID",  
                          to = "ENTREZID",  
                          multi.to = "first",  
                          multi.from = median 
                          ){  
  # ****************************************************** 
  # Function: RX_probe2gene
  #     
  #   This function is used to annotate gene names corresponding to probes 
  #   in an ExpressionSet or SummarizedExperiment object from a microarray
  #   study. 
  #   The assignment is performed using a preloaded annotation database.
  #
  # Input:
  # -----
  #   - Data (ExpressionSet or SummarizedExperiment):
  #         ExpressionSet or SummarizedExperiment object from the study.
  #   - anot_tb (data.frame, optional): 
  #         Annotation table. If provided, it is used instead of 
  #         the preloaded annotation database.
  #   - anotdb (AnnotationDbi, default "org.Hs.eg.db"): 
  #         Preloaded annotation package.
  #   - from (character, default "PROBEID"): 
  #         Current identifiers of the probes.
  #   - to (character, default "ENTREZID"): 
  #         Identifiers to which gene names will be assigned.
  #         
  #   - multi.to (character, default takes the first found value): 
  #         Handling of multiple matches with the target identifiers.
  #   - multi.from (function, default is the "median" function): 
  #         Handling of multiple matches with the source identifiers.
  # Output:
  # -----
  #   - ExpressionSet or SummarizedExperiment: 
  #         Modified object with gene names assigned to the probes.
  #
  # Dependencies:
  # -----
  #   library(AnnotationDbi)
  #
  # ******************************************************
  
  # Checking the type of object
  if (class(Data) == "SummarizedExperiment"){ 
    DataOut = Data
    exData = assay(Data)
  }else{
    DataOut = Data 
    exData = exprs(Data)
    #rowData(DataOut) = data.frame(rowData(DataOut))
  }
  # Annotation
  if(is.null(anot_tb)){
    Anot0 = AnnotationDbi::select(x = anotdb,
                                  keys = rownames(DataOut),
                                  columns = to,
                                  keytype = from)
  }else{
    Anot0 = anot_tb
  }
  
  # If there is more than one gene per probe, the first one is taken.
  ind1 = match(rownames(DataOut), Anot0[,from])
  Anot = Anot0[ind1,]
  
  # Eliminate probes that do not correspond to any gene. 
  Anot = Anot[which(! is.na(Anot[,to])),]
  
  # Join the annotation with the input data.
  Exp0 = data.frame(from = rownames(exData), exData)
  Exp1 = merge(Anot, Exp0, by.x = glue("{from}"), by.y = 1)
  
  # Take the name of the samples
  Samp = colnames(DataOut)
  
  # If there is more than one probe per gene use the ‘multi.from’ function.
  Exp2 = aggregate(Exp1[which(colnames(Exp1) %in% Samp)], 
                   by = list(X= Exp1[,to]),
                   multi.from)
  
  # We have in X the genes according to the nomenclature indicated.
  # Order as in the original.
  ord = match(unique(Exp1[,to]),  Exp2$X)
  Exp3 <- Exp2[ord,]
  ind2 = match(Exp3[,"X"], Anot[,to])
  Anot2 = Anot[ind2,] # Anot 2 is the final annotation.
  
  # Assign the values to the original matrix
  ind3 = match(rownames(DataOut), Anot2[,from]) 
  if (class(Data) == "SummarizedExperiment"){ 
    rowData(DataOut) = Anot2[ind3,]
    DataOut <- DataOut[which(!is.na(rowData(DataOut)$PROBEID),)]
    rownames(DataOut) = Exp3$X
    rownames(Exp3) = Exp3$X
    assay(DataOut, withDimnames=FALSE) <- as.matrix(Exp3[,-1])
  }else{
    fData(DataOut) = (Anot2[ind3,])
    DataOut <- DataOut[which(!is.na(fData(DataOut)$PROBEID),)]
    rownames(DataOut) = Exp3$X
    rownames(Exp3) = Exp3$X
    # We assign the values to the array
    exprs(DataOut) <- as.matrix(Exp3[,-1]) 
  }
  # Return the modified object
  return(DataOut)
}

#------------- RX_annot
RX_annot <- function(Data,  
                     db = org.Hs.eg.db,  
                     fromAnot="ENTREZID", 
                     toAnot = c("ENSEMBL","SYMBOL") 
                     ){
  # ****************************************************** 
  # Function: RX_annot
  #     
  #   This function is used to annotate identifiers in an 
  #   ExpressionSet or SummarizedExperiment object. It finds matches 
  #   between the current identifiers and the desired identifiers 
  #   using a preloaded annotation database.
  #
  # Input:
  # -----
  #   - Data (ExpressionSet or SummarizedExperiment): 
  #         Object containing the data to be annotated.
  #   - db (AnnotationDbPackage): 
  #         Preloaded annotation database package.
  #   - fromAnot (character, default is "ENTREZID"): 
  #         Current identifier of the data.
  #   - toAnot (character vector, default is "ENSEMBL" and "SYMBOL"): 
  #         Desired identifiers to be retrieved.
  #
  # Output:
  # -----
  #   - (ExpressionSet or SummarizedExperiment):
  #         Data object annotated with the desired identifiers.
  #
  # Dependencies:
  # -----
  #   library(AnnotationDbi)
  #   library(BiocGenerics)
  #   library(org.Hs.eg.db)
  #
  # ******************************************************
  
  Dataout = Data
  # Search for matches
  anot = AnnotationDbi::select(x = db,
                               keys = rownames(Dataout),
                               columns = toAnot,
                               keytype = fromAnot,
                               multiVals = "first")
  
  # Verify and manage multiple hits
  if(! isTRUE(all.equal(rownames(Dataout),anot[,fromAnot]))){  
    ind = BiocGenerics::match(rownames(Dataout), anot[,fromAnot]) 
    if (class(Data) == "SummarizedExperiment"){
      # Combinar con la información anterior
      rowData(Dataout) =cbind(anot[ind,],
                              rowData(Dataout)[which(!colnames(rowData(Dataout))
                              %in% c(fromAnot,toAnot))])
    }else{
      # Combine with the above information
      fData(Dataout) =cbind(anot[ind,],
                            fData(Dataout)[which(!colnames(fData(Dataout))
                            %in% c(fromAnot, toAnot))])
    }
  }
  return(Dataout)
}

# ~~~~~~~~~~~~ 03ExploratoryAnalysis ~~~~~~~~~~~~ #

#------------- RX_resumeBarPlot
RX_resumeBarPlot <- function(Data, 
                             Fac_1, 
                             Fac_2 = NULL, 
                             Fac_wrap = NULL, 
                             OneBar = TRUE, 
                             expandY = NULL, 
                             colors = NULL, 
                             alpha = c(0.2, 0.6),
                             ylab = "Number of samples",
                             xlab = Fac_1,
                             Fac_1_tit = Fac_1,
                             Fac_2_tit = Fac_2,
                             sep_width = 0.4,
                             Labels = TRUE,
                             numLabels = TRUE,
                             porcentLabels = TRUE,
                             numpos = 1.6,
                             porcentpos = 4,
                             size = 1,
                             text_size = 3.5,
                             porcent_size = 3
){
  # ****************************************************** 
  # Function: RX_resumeBarPlot
  #     
  #   This function is used to generate a bar plot summarizing the information 
  #   in an ExpressionSet or SummarizedExperiment object. The plot displays 
  #   the distribution of samples across groups defined by two factors
  #   (Fac_1 and Fac_2) and optionally a third factor (Fac_wrap).
  #
  # Input:
  # -----
  #   - Data (ExpressionSet or SummarizedExperiment): 
  #         Object containing the data to be visualized.
  #   - Fac_1 (character): 
  #         Name of the first factor grouping the samples on the x-axis 
  #         and defining the color of the bars.
  #   - Fac_2 (character, optional): 
  #         Name of the second factor that alters the transparency of the bars.
  #   - Fac_wrap (character, optional): 
  #         Name of the third factor that groups the data into subplots.
  #   - OneBar (logical, default TRUE): 
  #         Determines how samples are grouped.
  #         If TRUE, all levels of the second factor are shown within the same
  #         bar; if FALSE, each level of the factor has its own bar.
  #   - expandY (numeric, optional): 
  #         Expansion factor for the y-axis (e.g., 0.2 adds 20% more).
  #   - colors (vector, optional): 
  #         Dictionary of colors for the first factor.
  #   - alpha (numeric vector, default c(0.2, 0.6)): 
  #         Alpha values for the second factor 
  #         (bar transparency).
  #   - ylab (character, default "Number of samples"): 
  #         Label for the y-axis.
  #   - xlab (character, default Fac_1): 
  #         Label for the x-axis.
  #   - Fac_1_tit (character, default Fac_1): 
  #         Title of the legend for the first factor.
  #   - Fac_2_tit (character, default Fac_2): 
  #         Title of the legend for the second factor.
  #   - sep_width (numeric, default 0.4): 
  #         Spacing between bars.
  #   - Labels (logical, default TRUE): 
  #         Determines if labels are added to the bars.
  #   - numLabels (logical, default TRUE): 
  #         Adds labels with the number of samples per group.
  #   - porcentLabels (logical, default TRUE): 
  #         Adds labels with the percentage of samples per group 
  #         relative to the first factor.
  #   - numpos (numeric, default 1.6): 
  #         Position of the numerical labels on the y-axis relative to the top
  #         of the bar.
  #   - porcentpos (numeric, default 4): 
  #         Position of the percentage labels on the y-axis relative to the top
  #         of the bar.
  #   - size (numeric, default 1): 
  #         Border size of the bars.
  #   - text_size (numeric, default 3.5): 
  #         Size of the main text in the plot.
  #   - porcent_size (numeric, default 3.5): 
  #         Size of the percentage labels in the plot.
  # Output:
  # -----
  #   - (object of class "ggplot"):
  #         Bar plot generated with ggplot.
  #
  # Dependencies:
  # -----
  #   library(plyr)
  #   library(dplyr)
  #   library(ggplot2)
  #
  # ******************************************************
  
  # We convert the data
  Tabla <- Data %>%
    dplyr::group_by(RX_ifelse(is.null(Fac_1), NULL, Data[,Fac_1]),
                    RX_ifelse(is.null(Fac_2), NULL, Data[,Fac_2]),
                    RX_ifelse(is.null(Fac_wrap), NULL, Data[,Fac_wrap])) %>% 
    dplyr::summarise(Total = n())  
  
  colnames(Tabla) = c(RX_ifelse(is.null(Fac_1), NULL, "Fac_1"),
                      RX_ifelse(is.null(Fac_2), NULL, "Fac_2"),
                      RX_ifelse(is.null(Fac_wrap), NULL, "Fac_wrap"),
                      "Total")

  Tabla = Tabla %>%
    ddply(.variables = c("Fac_1"), transform, Porcentaje = round(Total/sum(Total)*100, 1)) %>%
    ddply(c("Fac_1", RX_ifelse(is.null(Fac_wrap), NULL, "Fac_wrap")), transform, Ypos = cumsum(Total)) 
  
  # Onebar
  if(isTRUE(OneBar)){
    # Plot
    BarPlot <- ggplot(Tabla, aes(x = RX_ifelse(is.null(Fac_1), "",Fac_1),
                                 y=Total,
                                 fill= RX_ifelse(is.null(Fac_1), "",Fac_1),
                                 color= RX_ifelse(is.null(Fac_1), "",Fac_1),
                                 alpha= RX_ifelse(is.null(Fac_2),
                                                  NULL, Fac_2))) + 
      geom_bar(width = sep_width,stat="identity", 
               position = position_stack(reverse = TRUE),
               size  = size 
      )  
    
    # Labels
    if(isTRUE(numLabels)){
      BarPlot <- BarPlot + 
        # Total
        geom_text(data = Tabla, aes( y = Ypos, label=glue("n = {Total}"), ), 
                  alpha = 1,
                  vjust= numpos,
                  color="black",
                  size=text_size) 
      if(isTRUE(porcentLabels)){
        BarPlot <- BarPlot +
          # Percentage
          geom_text(data = Tabla, aes( y = Ypos,
                                       label=glue("({Porcentaje} %)"), ), 
                    alpha = 1,
                    vjust=porcentpos,
                    color="black",
                    size=porcent_size) } 
    }
    # We resize the y-axis to 10 per cent more
    if(!is.null(expandY)){ 
      ymax0 = max(Tabla$Ypos)
      ymax = ymax0 +  ymax0 * expandY
      # we increase the Y-axis
      BarPlot <- BarPlot + ylim(c(0,ymax)) 
    }   
  }else{
    
    # Plot
    BarPlot <- ggplot(Tabla, aes(x = Fac_1,
                                 y=Total,
                                 fill= Fac_1,
                                 color= Fac_1,
                                 alpha= Fac_2)) + 
      geom_bar(width = sep_width,stat="identity", 
               position = position_dodge(), 
               size  = size
      )     
    # Labels
    if(isTRUE(Labels)){
      BarPlot <- BarPlot + 
        # Total
        geom_text(aes(label= glue("N = {Total}"), y = Total),
                  alpha = 1,vjust=1.6,                        
                  color="black", 
                  # hjust=0.5,
                  # Define text position and size
                  position = position_dodge2(sep_width),  
                  angle=0, 
                  size=text_size) + 
        # Percentage
        geom_text(aes(label= glue("({Percentage} %)"), y = Total),
                  alpha = 1,vjust=4,                        
                  color="black", 
                  #hjust=0.5,
                  # Define text position and size
                  position = position_dodge2(sep_width),  
                  angle=0, 
                  size=porcent_size) }
    # Resize the y-axis to 10 per cent larger
    if(!is.null(expandY)){ 
      ymax0 = max(Tabla$Total)
      ymax = ymax0 +  ymax0 * expandY
      # Increase Y-axis
      BarPlot <- BarPlot + ylim(c(0,ymax)) 
    }   
  }  
  # Group
  if(! is.null(Fac_wrap)){ 
    BarPlot <- BarPlot + facet_wrap(~Fac_wrap)}
  
  # Adjustment 
  # Colors
  if(! is.null(colors)){
    BarPlot <- BarPlot + 
      # We assign colours
      scale_fill_manual(values = colors) +  
      scale_color_manual(values = colors)
  }
  
  BarPlot <- BarPlot +  
    # Alpha
    scale_alpha_discrete(range = alpha) +
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill = "transparent"),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 0, size = 12)) +
    # Name of axles
    ylab(ylab) +
    xlab(xlab) + 
    # Legend
    guides(fill = guide_legend(title=Fac_1_tit,
                               order = 1),
           alpha = guide_legend(title=Fac_2_tit,
                                order = 2),
           color = "none")
  
  return(BarPlot)
}

#------------- RX_dataforPCA 
RX_dataforPCA <- function(data){
  # ****************************************************** 
  # Function: RX_dataforPCA
  #
  #   This function performs PCA computation on the provided data.
  #
  # Input:
  # -----
  #   - data (matrix or data frame):
  #         Data on which PCA is to be computed.
  #
  # Output:
  # -----
  #   - (object of class "prcomp"):
  #         Object containing the results of the PCA analysis.
  #
  # Dependencies:
  # -----
  #   library(stats)
  #
  # ******************************************************
  
  # Calculation of the transposed matrix
  TData = t(data)
  # PCA 
  DataPCA = prcomp(TData,scale.=TRUE)
  return(DataPCA)
}

#------------- RX_dataforPCA
RX_elimrows0 <- function(data){ 
  # ****************************************************** 
  # Function: RX_elimrows0
  #
  #   This function removes rows from a matrix or data frame that contain 
  #   only values equal to 0.
  #
  # Input:
  # -----
  #   - data (matrix or data frame):
  #         Data to be processed.
  #
  # Output:
  # -----
  #   - (matrix or data frame):
  #         Modified data (removing rows with only 0s).
  #
  # Dependencies:
  # -----
  #   none
  #
  # ******************************************************
  
  data = data[rowSums(data) > 0,]
  return(data)
}

#------------- RX_plotPCA
RX_plotPCA <- function(PCAData, phenoData, shade = 0.3,
                       shapes = c(16,17),
                       labels = NULL,
                       colors = NULL, 
                       Fac_wrap = NULL,
                       colorColumn = NULL,
                       shapeColumn = NULL, 
                       colorLegend = colorColumn, 
                       shapeLegend = shapeColumn,
                       alpha = 0.8,
                       size = 5
                       ){ 
  # ****************************************************** 
  # Function: RX_plotPCA
  #
  #   This function creates a PCA plot from a PCA object and phenotypic data.
  #
  # Input:
  # -----
  #   - PCAData (PCA object):
  #         PCA object obtained from the prcomp function.
  #   - phenoData (data frame):
  #         Phenotypic data used to add information to the plot points.
  #   - shade (numeric, default 0.3):
  #         Shading value for background points.
  #   - shapes (numeric vector, default c(16,17)):
  #         Vector of shape codes for the plot points.
  #   - labels (character, optional):
  #         Name of the variable to be used for point labels.
  #   - colors (vector of colors, optional):
  #         Vector of colors for the plot points.
  #   - Fac_wrap (character, optional):
  #         Name of the variable used to separate points into multiple panels.
  #   - colorColumn (character, optional):
  #         Name of the column used to assign colors to the points.
  #   - shapeColumn (character, optional):
  #         Name of the column used to assign shapes to the points.
  #   - colorLegend (character, default "colorColumn"):
  #         Name of the column used for the color legend title. If NULL, the color legend is not shown.
  #   - shapeLegend (character, default "colorColumn"):
  #         Name of the column used for the shape legend title. If NULL, the shape legend is not shown.
  #   - alpha (numeric, default 0.8):
  #         Transparency value for the plot points.
  #   - size (numeric, default 5):
  #         Size of the plot points.
  #
  # Output:
  # -----
  #   - (object of class "ggplot"):
  #         PCA plot generated with ggplot.
  #
  # Dependencies:
  # -----
  #   - library(ggplot2)
  #   - library(ggrepel)
  #   - library(stats)
  #
  # ******************************************************
  
  # Data conversion 
  inputData <- data.frame(PCAData$x[,c(1,2)], phenoData) 
  
  # Creation of the ggplot object
  PCAplot = ggplot(inputData) 
  
  # Factorization
  if (! is.null(Fac_wrap)){
    PCAplot = PCAplot +
      # Background points
      geom_point(data= inputData[,which(colnames(inputData)!=Fac_wrap)],
                 aes(PC1, PC2,
                     colour  = "grey",
                     shape = Sex,
                     alpha = shade,
                     size = 5))
    } 
  # Adding PCA points    
  PCAplot = PCAplot +
    geom_point(aes(PC1, PC2,
                     colour  = RX_ifelse(is.null(colorColumn),
                                         NULL,inputData[,colorColumn]),
                     shape = RX_ifelse(is.null(shapeColumn),
                                       NULL,inputData[,shapeColumn])
                   ), 
               alpha = alpha,
               size = size)  + 
      # Theme
      theme(plot.background = element_rect(fill="transparent"),
            panel.background = element_rect(fill = "transparent"),
            axis.line = element_line(colour = "black")) +
      # Colors
      scale_color_manual(values = colors) +
      # Shape
      scale_shape_manual(values = shapes) 
    
    if(! is.null(Fac_wrap)){
      # Tissue division
      PCAplot = PCAplot +
        facet_wrap(facets = c(Fac_wrap), ncol = 2, scales = "free")} 
  
      # Title legends
      # Legend
  PCAplot = PCAplot + guides(color = RX_ifelse((is.null(colorLegend) || 
                                                  colorLegend == "none"),
                                               "none",
                               guide_legend(title = colorLegend,
                                            order = 1,
                                            override.aes=list(size = 4))),
             shape = RX_ifelse((is.null(shapeLegend) || 
                                  shapeLegend == "none"), "none",
                               guide_legend(title = shapeLegend,
                                            order = 3,
                                  override.aes=list(size = 4))),
             alpha = "none",
             size = "none"
             )  + 
      # Titulos de los ejes
      xlab(glue("PCA1 ({round(summary(PCAData)$importance[2,1], 4) *100}% of explained variance)"))+
      ylab(glue("PCA2 ({round(summary(PCAData)$importance[2,2], 4) *100}% of explained variance)"))
  
    
  # Scale
  escalx = max(inputData$PC1) - min(inputData$PC1)
  escaly = max(inputData$PC2) - min(inputData$PC2) 
  # Points labels
  if(!is.null(labels)){
    PCAplot <- PCAplot + 
      geom_text_repel(aes(x = PC1, y=PC2, label = phenoData$Sample),
                      arrow = arrow(length = unit(0.01, "npc"), 
                                    type = "closed", ends = "last"), 
                      segment.size  = 0.2,
                      segment.color = "grey50",
                      direction     = "x", 
                      max.overlaps = 25,
                      size= 2.5,
                      nudge_y = 0.02 * escaly,
                      nudge_x = 0.025 * escalx
                      )  
  }
  return(PCAplot)
}

#------------- RX_multiplotPCA
RX_multiplotPCA <- function(Data,
                            multiplot = FALSE,
                            factor = NULL,
                            ...){  
  # ****************************************************** 
  # Function: RX_multiplotPCA
  #
  #   This function removes rows from a matrix or data frame that contain
  #   only zero values.
  #
  # Input:
  # -----
  #   - Data (ExpressionSet or SummarizedExperiment):
  #         Expression data for performing the PCA analysis.
  #   - multiplot (logical, default FALSE):
  #         Indicates whether to create a separate PCA for each level of the factor
  #         or a single PCA for all the data.
  #   - factor (character, optional):
  #         Name of the column used to separate the data into multiple PCAs.
  #         Used only if multiplot = TRUE.
  #   - ...:
  #         Other arguments passed to the RX_plotPCA function.
  #
  # Output:
  # -----
  #   - (object of class "ggplot" or list of "ggplot" objects):
  #         PCA plot(s) generated with ggplot.
  #
  # Dependencies:
  # -----
  #   - library(glue)
  #   - library(ggplot2)
  #   - library(ggrepel)
  #   - library(stats)
  #
  # ******************************************************
  
  if(isFALSE(multiplot)){
    # We create a single PCA 
    PCAData = 
      RX_dataforPCA(
        RX_elimrows0(
          RX_ifelse(class(Data) == "SummarizedExperiment",
                    assay(Data),
                    exprs(Data)
          )
        )
      )
    
    phenoData = 
      RX_ifelse(class(Data) == "SummarizedExperiment",
                colData(Data),
                pData(Data)
      ) 
    # Plot
    PCAplot = RX_plotPCA(PCAData, phenoData, Fac_wrap = factor, ...)
    return(PCAplot)
  }else{
    # We create a PCA for each level of the factor
    # Levels to be taken
    Fac = glue("Data${factor}")
    Levs = eval(parse(text=Fac))
    # Division
    T1 = "Data[,which(Levs == x)]" 
    
    # Martix
    PCADatas =
      sapply(levels(Levs), function(x) 
        RX_dataforPCA(
          RX_elimrows0(
            RX_ifelse(class(eval(parse(text = T1))) == "SummarizedExperiment",
                      assay(eval(parse(text = T1))),
                      exprs(eval(parse(text = T1)))
            )
          )
        ), simplify = FALSE, USE.NAMES = TRUE
      ) 
    # Pheno data
    phenoDatas =
      sapply(levels(Levs), function(x) 
        data.frame(
          RX_ifelse(class(eval(parse(text = T1))) == "SummarizedExperiment",
                    colData(eval(parse(text = T1))),
                    pData(eval(parse(text = T1)))
          )
        ), simplify = FALSE, USE.NAMES = TRUE
      ) 
    # Plots
    PCAplots = 
      lapply(levels(Levs), function(x)
        RX_plotPCA(PCADatas[[x]], phenoDatas[[x]], Fac_wrap = factor, ...)
      )
    'PCAplotsOut = do.call("ggarrange", c(PCAplots,
                                         common.legend = TRUE, 
                                         legend = "right"))'
    return(PCAplots)
  }
}

#------------- RX_dataforboxplot
RX_dataforboxplot <- function(data, 
                              transformacion = NULL, 
                              phenoData = NULL,
                              ncol=1 ){
  # ****************************************************** 
  # Function: RX_dataforboxplot
  #
  #   This function removes rows from a matrix or data frame that contain
  #   only zero values.
  #
  # Input:
  # -----
  #   - data (matrix):
  #         Expression matrix.
  #   - transformation (character, optional):
  #         Type of transformation to apply to the data.
  #         Can be "log2" or "log10". Default is NULL, meaning no transformation will be applied.
  #   - phenoData (data frame or matrix, optional):
  #         Phenotypic data used to add additional information to the boxplots.
  #   - ncol (integer, default 1):
  #         Index of the column in phenoData that contains the sample names.
  #
  # Output:
  # -----
  #   - (object of class "data.frame"):
  #         Data prepared for creating a boxplot.
  #
  # Dependencies:
  # -----
  #   - library(dplyr)
  #   - library(tidyr)
  #
  # ******************************************************

  # Apply transformation if specified
  if(! is.null(transformacion) ){
    if(transformacion == "log2"){
      data = log2(data)
    }else{
      if(transformacion == "log10"){
        data = log10(data)
      }
    }}
  # Convert data to data.frame
  data = data.frame(data)
  
  # Prepare the data for the pivot_longer function
  Dt = data %>% pivot_longer(cols = c(colnames(data)),
                             names_to = "Sample",
                             values_to = "exp")
  
  # Sort data according to originals
  Dt$Sample = factor(Dt$Sample,
                     levels = colnames(data))
  Dt = Dt[(order(Dt$Sample)),]
  
  # Join phenotypic data if specified
  if(! is.null(phenoData)){
    Dt = merge(x = Dt, y= phenoData, 
               by.x = "Sample", by.y = ncol,
               sort = FALSE)}
  
  return(Dt)
}

#------------- RX_BoxPlot 
RX_BoxPlot <- function(Data,
                       Fac_fill = "Group",
                       Fac_alpha = "Sex",
                       Fac_wrap= "Tissue",
                       colors = c("#90C432", "#2494B5"),
                       alpha = c(0.95, 0.45),
                       shape = 19,
                       dir = "h",# 8
                       ncol = 1
){
  # ****************************************************** 
  # Function: RX_BoxPlot
  #
  # Creates a boxplot.
  #
  # Input:
  # -----
  #   - Data (data frame):
  #         Data used to create the boxplot.
  #   - Fac_fill (character, default "Group"):
  #         Name of the column used for filling the boxplots.
  #   - Fac_alpha (character, default "Sex"):
  #         Name of the column used for the transparency of the boxplots.
  #   - Fac_wrap (character, default "Tissue"):
  #         Name of the column used to split the boxplots into multiple panels.
  #   - colors (vector of colors, default c("#90C432", "#2494B5")):
  #         Colors used for filling the boxplots.
  #   - alpha (numeric vector, default c(0.95, 0.45)):
  #         Transparency values used for the boxplots.
  #   - shape (integer, default 19):
  #         Shape code used for outlier points in the boxplots.
  #   - dir (character, default "h"):
  #         Direction for splitting the panels. "h" for horizontal and
  #         "v" for vertical.
  #   - ncol (integer, default 1):
  #         Number of columns used for arranging the panels.
  #
  # Output:
  # -----
  #   - (object of class "ggplot"):
  #         Boxplot created with ggplot.
  #
  # Dependencies:
  # -----
  #   - library(ggplot2)
  #
  # ******************************************************
  
  BoxPlot = (
    # ggplot object
    ggplot(Data,
           aes(x = Sample, y = exp)) +
      # Boxplot
      geom_boxplot(aes(fill = Data[,Fac_fill],
                       alpha = Data[,Fac_alpha]),
                   outlier.colour = "black", outlier.shape = shape,
                   outlier.alpha = 1)+ 
      # Theme
      theme(plot.background = element_rect(fill="transparent"),
            panel.background = element_rect(fill = "transparent"),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 1, size = 8)) +
      # Fill colours
      scale_fill_manual(values = colors) +  
      # Intensity of colours
      scale_alpha_discrete(range = alpha) + 
      # Tissue division
      facet_wrap(~get(Fac_wrap), ncol = ncol, scales = "free", dir = dir) +
      # Title legends
      # Legend
      guides(fill = guide_legend(title="Grupo",
                                 order = 1),
             alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),100,0,
                                                           alpha=alpha),
                                                  colour="black"),
                                title = "Sexo",
                                order = 2)) +
      # Title axles
      xlab("Samples") +
      ylab("log2(expression)") +
      # Labels 2 levels x-axis
      scale_x_discrete(guide = guide_axis(n.dodge = 2))
  )
  return(BoxPlot)
}


#------------- RX_Clustering 
RX_Clustering <- function(expressData,
                          method = c("euclidean", "correlation")
                          ){
  # ****************************************************** 
  # Function: RX_Clustering
  #
  #   This function computes hierarchical clustering using either correlation or Euclidean distance.
  #
  # Input:
  # -----
  #   - expressData (matrix or data frame):
  #         Expression data used for clustering.
  #   - method (character, default c("euclidean", "correlation")):
  #         Method used to calculate the distance.
  #         It can be "euclidean" for Euclidean distance or
  #         "correlation" for Pearson's correlation coefficient.
  #
  # Output:
  # -----
  #   - (object of class "hclust"):
  #         Hierarchical clustering object obtained.
  #
  # Dependencies:
  # -----
  #   - library(stats)
  #
  # ******************************************************
  
  # Herarchical clustering and data transformation
  if(method == "euclidean"){
    Tans <- t(expressData)
    Dist   <- dist(Tans, method = "euclidean")
    hc  <- hclust(Dist)  
  }else{
    Correlacion <- cor(expressData)
    Dist <- as.dist((1 - Correlacion) / 2)
    hc <- hclust(Dist, method = "complete") 
  }
  return(hc)
}

#------------- RX_datatodendo 
RX_datatodendo <- function(hc, 
                           phenoData = NULL, 
                           ncol = 1){ # Add phenotypic information to a dendogram
  
  # ****************************************************** 
  # Function: RX_datatodendo
  #
  #   This function adds phenotypic information to a dendrogram.
  #
  # Input:
  # -----
  #   - hc (object of class "hclust"):
  #         Hierarchical clustering object obtained previously.
  #   - phenoData (data frame, optional):
  #         Phenotypic data used to add information to the dendrogram nodes.
  #   - ncol (numeric, default 1):
  #         Column index in phenoData that matches the node names 
  #         in the hc object.
  #
  # Output:
  # -----
  #   - (object of class "dendro"):
  #         Dendrogram object with added phenotypic information.
  #
  # Dependencies:
  # -----
  #   - library(ggdendro)
  #
  # ******************************************************
  
  
  # We convert the hierarchical clustering object to a dendrogram.
  dendoData0 <- as.dendrogram(hc)
  dendoData <- dendro_data(dendoData0)
  
  # We add the phenotypic data
  if (!is.null(phenoData)){
    dendoData$labels <- data.frame(merge(dendoData$labels, phenoData, 
                                         by.x = "label", by.y=ncol,
                                         sort = FALSE))
  }
  # We return the object
  return(dendoData)
}

#------------- RX_dendoPlot 
RX_dendoPlot <- function(dendoData,
                         dir = c("lr", "rl", "tb", "bt"),
                         colorPal = NULL,
                         colorOutPal = NULL,
                         colorOutColumn = NULL,
                         branchSize = NULL,
                         labelSize = 2,
                         moveLabel = 0.001,
                         expandY = 0.1,
                         expandYlbs = 0,
                         colorColumn = NULL,
                         leaves = FALSE,
                         shapeColumn = NULL,
                         leavesSize = 2,
                         leavesStroke = 1,
                         leavesShapes = NULL,
                         colorLegend = colorColumn,
                         colorOutLegend = colorOutColumn,
                         shapeLegend = shapeColumn
                         
){
  # ****************************************************** 
  # Function: RX_dendoPlot
  #
  #   This function draws a dendrogram from the data provided in dendoData.
  #
  # Input:
  # -----
  #   - dendoData (object of class "dendro"):
  #         Dendrogram object with the necessary information to draw the dendrogram.
  #   - dir (character, default "lr"):
  #         Direction of the dendrogram drawing. It can be "lr" (left to right),
  #         "rl" (right to left), "tb" (top to bottom), or "bt" (bottom to top).
  #   - colorPal (vector of characters, optional):
  #         Color palette used for the points on the leaves of the dendrogram.
  #   - colorOutPal (vector of characters, optional):
  #         Color palette used for the branches and segments of the dendrogram.
  #   - colorOutColumn (character, optional):
  #         Column name in dendoData containing the variable used for the color of the branches and segments.
  #   - branchSize (numeric, optional):
  #         Size of the branches of the dendrogram.
  #   - labelSize (numeric, default 2):
  #         Size of the labels in the dendrogram.
  #   - moveLabel (numeric, default 0.001):
  #         Vertical shift of the labels in the dendrogram.
  #   - expandY (numeric, default 0.1):
  #         Vertical expansion factor of the dendrogram.
  #   - expandYlbs (numeric, default 0):
  #         Vertical expansion factor of the labels relative to the total length of the dendrogram.
  #   - colorColumn (character, optional):
  #         Column name in dendoData containing the variable used for the color of the points on the leaves.
  #   - leaves (logical, default FALSE):
  #         Indicates whether to show the points on the leaves of the dendrogram.
  #   - shapeColumn (character, optional):
  #         Column name in dendoData containing the variable used for the shape of the points on the leaves.
  #   - leavesSize (numeric, default 2):
  #         Size of the points on the leaves of the dendrogram.
  #   - leavesStroke (numeric, default 1):
  #         Stroke width of the points on the leaves of the dendrogram.
  #   - leavesShapes (vector of characters, optional):
  #         Shapes used for the points on the leaves of the dendrogram.
  #   - colorLegend (character, optional):
  #         Column name in dendoData containing the variable used for the color legend.
  #   - colorOutLegend (character, optional):
  #         Column name in dendoData containing the variable used for the color legend of the branches.
  #   - shapeLegend (character, optional):
  #         Column name in dendoData containing the variable used for the shape legend.
  #
  # Output:
  # -----
  #   - (object of class "ggplot"):
  #         ggplot object representing the dendrogram.
  #
  # Dependencies:
  # -----
  #   - library(ggnewscale)
  #   - library(ggplot2)
  #
  # ******************************************************
  
  # Parameters
  dir <- match.arg(dir)
  ybreaks   <- pretty(segment(dendoData)$y, n = 5)
  ymax      <- max(segment(dendoData)$y)
  
  # Basic Dendrogram
  DendoPlot <- ggplot(segment(dendoData)) + 
    geom_segment(aes(x = x, y = y, 
                     xend = xend,
                     yend = yend, 
                     #lineend="round",
                     size =  branchSize)) + 
    if(! is.null(colorColumn)){
      geom_segment(data = dendoData$segments %>%
                     filter(yend == 0) %>%
                     left_join(dendoData$labels, by = "x"), 
                   aes(x=x, y=y.x, xend=xend, yend=yend, 
                       color = dendoData$labels[,colorColumn]))
    }
  #coord_flip() + 
  #scale_y_reverse(expand = c(0.2, 0))
  
  DendoPlot <- DendoPlot + scale_x_continuous(breaks = NULL)
  if (dir %in% c("rl", "lr")) {
    DendoPlot <- DendoPlot + coord_flip()
  }
  if (dir %in% c("bt", "lr")) {
    DendoPlot <- DendoPlot + scale_y_reverse(breaks = ybreaks)
  } else {
    DendoPlot <- DendoPlot + scale_y_continuous(breaks = ybreaks)
    moveLabel <- -(moveLabel)
  }
  
  # Labels
  # Se params
  #{
  nLabs = nrow(dendoData$labels)
  angle <- rep(0, nLabs)
  hjust <- 0
  if (dir %in% c("tb", "bt")) {angle <- angle + 45 }
  if (dir %in% c("tb", "rl")) {hjust <- 1 }
  
  labelParams = list(angle = angle, hjust = hjust, vjust = 0.5)
  dendoData$labels$angle <- labelParams$angle 
  #}
  
  # Add labels
  DendoPlot <- DendoPlot + 
    geom_text(data = data.frame(label(dendoData)),
              aes(x = x,
                  y = (y + (ymax * expandYlbs)),
                  label = label,
                  colour = RX_ifelse(is.null(colorColumn), NULL,
                                     label(dendoData)[,colorColumn]),
                  angle = angle),
              vjust = labelParams$vjust,
              hjust = labelParams$hjust,
              nudge_y = ymax * moveLabel,
              size = labelSize,
              show.legend = FALSE)
  # Assign the colour palette
  DendoPlot <- DendoPlot + guides(color = "none")
  
  if (!is.null(colorPal)) {
    DendoPlot <- DendoPlot + scale_color_manual(values = colorOutPal)
  } 
  
  # We add the points  
  if (! isFALSE(leaves)){
    
    DendoPlot <- DendoPlot +
      new_scale_color() +
      geom_point(data = dendoData$labels,  
                 aes(x = seq(dendoData$labels$label),
                     y    = 0.00,
                     shape = RX_ifelse(is.null(shapeColumn), 
                                       NULL,label(dendoData)[,shapeColumn]),
                     fill = RX_ifelse(is.null(colorColumn), 
                                      NULL,label(dendoData)[,colorColumn]), 
                     color = RX_ifelse(is.null(colorOutColumn), 
                                       NULL,label(dendoData)[,colorOutColumn]),
                     alpha = NULL), 
                 size = leavesSize,
                 stroke = leavesStroke, 
                 show.legend = TRUE,
      ) 
    
    if(!is.null(leavesShapes)){
      DendoPlot <- DendoPlot + scale_shape_manual(values = leavesShapes)
    }
    if (!is.null(colorOutPal)) {
      DendoPlot <- DendoPlot + scale_color_manual(values = colorOutPal)
    }
    
    if (!is.null(colorPal)) {
      DendoPlot <- DendoPlot + scale_fill_manual(values = colorPal)
    }}
  
  # We modify the length
  ylim <- ymax * expandY#-round(ymax * expandY, 1)
  DendoPlot <- DendoPlot + expand_limits(y = ylim)
  
  # Legend
  
  DendoPlot <- DendoPlot + 
    guides(color = RX_ifelse(colorOutLegend == "none", "none",
                             guide_legend(title = colorOutLegend,
                                          override.aes=list(colour=colorOutPal),
                                          order = 1)),
           fill = RX_ifelse(colorLegend == "none", "none",
                            guide_legend(title = colorLegend,
                                         override.aes=list(colour=colorPal),
                                         order = 2)),
           shape = RX_ifelse(shapeLegend == "none", "none",
                             guide_legend(title = shapeLegend,
                                          order = 3))
  ) + 
    theme_dendro() + 
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill = "transparent"),
          axis.line = element_line(colour = "transparent"))
  
  return(DendoPlot)
}

#------------- RX_multiplotDendo 
RX_multiplotDendo <- function(Data, 
                              method = "correlation", 
                              multiplot = FALSE, 
                              factor = NULL,
                              wrap = FALSE, 
                              ...){ 
  # ****************************************************** 
  # Function: RX_multiplotDendo
  #
  #   This function creates dendrograms from the provided data and organizes
  #   them into one or multiple figures, depending on the specified parameters.
  #
  # Input:
  # -----
  #   - Data (object of class "SummarizedExperiment" or "matrix"):
  #         Data used to build the dendrograms.
  #   - method (character, default "correlation"):
  #         Method used to calculate the distance matrix for hierarchical clustering. 
  #         Can be "correlation", "euclidean", "manhattan", or other methods available in the 
  #         "dist" function.
  #   - multiplot (logical, default FALSE):
  #         Indicates whether to create individual dendrograms for each level of the specified factor.
  #   - factor (character, optional):
  #         Name of the factor used to create individual dendrograms. Only used if "multiplot" is TRUE.
  #   - wrap (logical, default FALSE):
  #         Indicates whether the individual dendrograms should be wrapped in a grid layout.
  #   - ... :
  #         Other arguments passed to the "RX_dendoPlot" function to customize the dendrograms.
  #
  # Output:
  # -----
  #   - (object of class "ggplot"):
  #         ggplot object representing the dendrogram.
  #
  # Dependencies:
  # -----
  #   - library(glue)
  #   - library(ggnewscale)
  #   - library(ggplot2)
  #   - library(stats)
  #
  # ******************************************************
  
  if(isFALSE(multiplot)){
    # We create a single Dendogram  
    # Data
    expressData = RX_ifelse(class(Data) == "SummarizedExperiment",
                            assay(Data),
                            exprs(Data)
    )
    phenoData = RX_ifelse(class(Data) == "SummarizedExperiment",
                          colData(Data),
                          pData(Data)
    ) 
    hc = RX_Clustering(expressData, method)
    dendoData = RX_datatodendo(hc, phenoData, ncol = "Sample")  
    # Plot
    Dendoplot = RX_dendoPlot(dendoData, ...) 
    return(list(Data = dendoData, Plot = Dendoplot))
    
  }else{
    # We create one for each level of the factor
    # Levels
    Fac = glue("Data${factor}")
    Levs = eval(parse(text=Fac))
    # Division
    T1 = "Data[,which(Levs == x)]" 
    
    # Martix
    DendoDatas =
      sapply(levels(Levs), function(x) 
        RX_datatodendo(
          RX_Clustering(
            RX_ifelse(class(eval(parse(text = T1))) == "SummarizedExperiment",
                      assay(eval(parse(text = T1))),
                      exprs(eval(parse(text = T1)))
            ), method
          ), 
          data.frame(
            RX_ifelse(class(eval(parse(text = T1))) == "SummarizedExperiment",
                      colData(eval(parse(text = T1))),
                      pData(eval(parse(text = T1)))
            )
          ), ncol = "Sample"
          
        ), simplify = FALSE, USE.NAMES = TRUE
      )
    # Plots
    Dendoplots = 
      lapply(levels(Levs), function(x)
        RX_dendoPlot(DendoDatas[[x]], ...) +
          if(wrap == TRUE){
            # Tissue division
            facet_wrap(~Tissue, nrow = 3, scales = "free") 
          }else{
            NULL
          }
      )
    return(list(Data = DendoDatas, Plot = Dendoplots))
  }
}

#------------- RX_DendoRelations 
RX_DendoRelations <- function(DendoData1,
                              DendoData2,
                              size= c(0.65,10,0)){
  # ****************************************************** 
  # Function: RX_DendoRelations
  #
  #   Creates a visual representation of the relationships between two dendrograms.
  #
  # Input:
  # -----
  #   - DendoData1 (object of class "dendro"):
  #         Data for the first dendrogram.
  #   - DendoData2 (object of class "dendro"):
  #         Data for the second dendrogram.
  #   - size (numeric vector, default c(0.65, 10, 0)):
  #         The relative size of the different parts of the plot. The vector
  #         should have three elements representing the size of the top panel,
  #         the size of the line panel, and the size of the bottom panel, respectively.
  #
  # Output:
  # -----
  #   - (object of class "ggplot"):
  #         ggplot object representing the dendrogram relationship.
  #
  # Dependencies:
  # -----
  #   - library(ggplot2)
  #   - library(ggpubr)
  #
  # ******************************************************
  
  TH = theme(plot.background = element_rect(fill="transparent"),
             panel.background = element_rect(fill = "transparent", 
                                             colour = "transparent"))
  
  idx <- data.frame(y1 = 1:nrow(DendoData1$labels),
                    y2 = match(DendoData1$labels$label, DendoData2$labels$label)) 
  plt <- ggplot() +
    geom_segment(data     = idx, 
                 aes(x    = 0,
                     y    = y1,
                     xend = 1,
                     yend = y2),
                 color    = "black") + 
    theme_dendro() + TH
  plt <- ggarrange(ggplot() + theme_void() ,plt, 
                   ggplot()+ theme_void() , ncol=1, heights = size) +
    guides(fill  ="none", shape="none")
  return(plt)
}


# ~~~~~~~~~~~~ 04DifferentialExpression ~~~~~~~~~~~~ #

#------------- RX_contrast
RX_contrast <- function(expmatrix,
                        design,
                        C,
                        seed=1808){
  # ****************************************************** 
  # Function: RX_contrast
  #
  #   This function performs a differential expression analysis using 
  #   the "limma" package. 
  #   It computes the results for a given expression dataset "expmatrix", 
  #   the design matrix "design", and the specific contrast "C".
  #   Finally, it adjusts the p-values using the "BH" method and 
  #   returns the results.
  #
  # Input:
  # -----
  #   - expmatrix (data.frame or matrix):
  #         The expression dataset used in the contrast analysis.
  #   - design (matrix):
  #         The experimental design matrix used in the contrast analysis.
  #   - C (character):
  #         The specific contrast to compute.
  #   - seed (numeric value, default 1808):
  #         The seed used for generating random numbers in the analysis.
  #
  # Output:
  # -----
  #   - (object of class "MArrayLM"):
  #       Results from the contrast analysis.
  #
  # Dependencies:
  # -----
  #   - library(limma)
  #
  # ******************************************************
  
  set.seed(seed)
  # Contrast matrix
  contMatrix <- makeContrasts(contrasts = C, levels = design)
  # Fit a linear model
  fit <- lmFit(expmatrix, design)
  fit2 <- contrasts.fit(fit, contMatrix)
  # Differential expression analysis
  Res <- eBayes(fit2)
  Res$TopTab <- topTable(Res, number = Inf, adjust.method = "BH", sort.by = "none") #, sort.by = "logFC"
  
  return(Res)
}

#------------- RX_DiffExp
RX_DiffExp <- function(phenoData, 
                       expressData, 
                       C, 
                       var="Group",
                       studyType = "Array",
                       covar = NULL,
                       seed=1808){
  # ****************************************************** 
  # Function: RX_DiffExp
  #
  #   This function performs a differential expression analysis using 
  #   the "limma" package. 
  #   It takes as input a `phenoData` object containing phenotypic information 
  #   of the samples, and an `expressData` object containing the expression data.
  #   It computes the contrasts specified in "C" for the expression data 
  #   "expressData" and the phenotypic data "phenoData".
  #   The analysis can be performed for "Array" or "RNA-seq" type studies, 
  #   and covariates can be included in the design.
  #   It returns a list of objects containing the results for each contrast.
  #
  # Input:
  # -----
  #   - phenoData (data.frame):
  #         The phenotypic data.
  #   - expressData (data.frame or matrix):
  #         The expression data.
  #   - C (character vector):
  #         Vector of contrasts to compute.
  #   - var (character or vector of characters, default "Group"):
  #         The name or names of the variables in "phenoData" to be used 
  #         to construct the contrasts.
  #         If two names are provided, an interaction between the variables 
  #         will be created.
  #   - studyType (character, default "Array"):
  #         The study type, which can be "Array" or "RNA-seq". If "RNA-seq", 
  #         the "voom" method will be applied to the expression data.
  #   - covar (character, optional):
  #         The name of the variable in "phenoData" to be used as a covariate 
  #         in the design.
  #   - seed (numeric value, default 1808):
  #         The seed used for generating random numbers in the analysis.
  #
  # Output:
  # -----
  #   - (list)
  #     List of results by contrast. Each object is of class "MArrayLM" 
  #     and contains the differential expression analysis results for a 
  #     specific contrast.
  #
  # Dependencies:
  # -----
  #   - library(limma)
  #   - library(Biobase)
  #   - library(stringr)
  #
  # ******************************************************
  
  set.seed(seed)
  
  if(length(var) == 2){
    # Create an interaction contrast between the two phenoData variables
    contrast = interaction(phenoData[,var[1]], phenoData[,var[2]])
  }else{
    contrast = phenoData[,var]
  }
  # Evaluating the presence of covariates
  if(is.null(covar) | ! is.factor(phenoData[, covar])){
    # Create a design matrix without covariates
    design = model.matrix(~0 + contrast)
    colnames(design) <- str_replace(levels(contrast) , "-", "_")
  }else{
    # Create a design matrix with the contrast and the covariate
    cov = factor(make.names(phenoData[, covar]))
    design = model.matrix(~ 0 + contrast + cov)
    colnames(design) <- c(str_replace(levels(contrast) , "-", "_"),
                          levels(cov)[-2])
  }
  # If studyType is ‘RNA-seq’, apply the voom transformation.
  if (studyType == "RNA-seq"){
    #library(edgeR)
    #d0 <- DGEList(expressData)
    #d0 <- calcNormFactors(d0, method="TMM")
    expressData = voom(expressData, design = design, plot = F)
  }
  # Perform differential expression analysis for each contrast.
  Res <- sapply(C, function(x) RX_contrast(expressData, design, x),
                simplify = FALSE, USE.NAMES = TRUE)
  return(Res)
}


#------------- RX_DiffExpFinal
RX_DiffExpFinal <- function(Data,
                            var1="Group",
                            var2="Sex",
                            seed = 1808,
                            ...){
  # ******************************************************** 
  # Function: RX_DiffExpFinal
  #
  #   This function performs a set of differential expression analyses for the 
  #   scenarios of a study using the "limma" package. 
  #   It performs differential expression analysis for different tissues in a 
  #   dataset. For each of these, it calculates the contrasts specified with 
  #   the "RX_DiffExp" function and combines the results into a list. 
  #
  # Input:
  # -----
  #   - Data (ExpressionSet or SummarizedExperiment object):
  #         The dataset containing both the expression data and the phenotype data.
  #   - var1 (character, default "Group"):
  #         The name of the variable in the phenotype data used to construct the contrasts.
  #   - var2 (character, default "Sex"):
  #         The name of the variable in the phenotype data used to construct the contrasts.
  #   - seed (integer, default 1808):
  #         The seed used for generating random numbers.
  #   - ... :
  #         Other arguments that will be passed to the "RX_DiffExp" function.
  #
  # Output:
  # -----
  #   - List of differential expression results.
  #     For each tissue, study, and contrast, we obtain an object of class MArrayLM.
  #
  # Dependencies:
  # -----
  #   - library(SummarizedExperiment)
  #   - library(Biobase)
  #   - library(glue)
  #
  # ********************************************************

  set.seed(seed)
  # We determine the provenance of the expression data
  if(class(Data) == "ExpressionSet"){
    phenoData = pData(Data)
    expressData = exprs(Data)
    studyType = "Array"
  }else{
    phenoData = colData(Data)
    studyType = "RNA-seq"
    expressData = assay(Data)      
  }
  # We take the variables
  Group = str_replace(levels(phenoData[,var1]) , "-", "_")
  Sex = levels(phenoData[,var2])
  Tissue = levels(phenoData$Tissue)
  # We separate by tissue
  phenoDatas = sapply(Tissue, function(x)
    phenoData[which(phenoData$Tissue == x),],
                      simplify = FALSE, USE.NAMES = TRUE)
  expressDatas = sapply(Tissue, function(x)
    expressData[,which(phenoData$Tissue == x)],
                        simplify = FALSE, USE.NAMES = TRUE)
  
  V1 = str_replace(levels(phenoData[,var1]) , "-", "_")
  V2 = levels(phenoData[,var2]) 
  # Contrasts
  C1 = combn(rev(Group), m=2, FUN = paste,collapse=" - ", simplify = TRUE)
  
  C2 = sapply(Sex, 
              function(x) 
                combn(levels(interaction(Group, x))[match(rev(Group), strsplit2(levels(interaction(Group, x)),
                                                                                          split = ".", fixed = TRUE)[,1])],
                                m=2, FUN = paste,collapse=" - ",
                                simplify = TRUE), simplify = F)
  
  C3 = sapply(seq(C1), function(i) glue("({C2$M[i]}) - ({C2$F[i]})"))
  C3 = c(unlist(C2), C3, use.names = FALSE)
  
  
  Tt = sapply(Tissue, function(x) c(RX_DiffExp(phenoDatas[[x]], 
                                               expressDatas[[x]] , 
                                               C1, var=var1, 
                                               studyType = studyType, ...),
                                    RX_DiffExp(phenoDatas[[x]], 
                                               expressDatas[[x]], 
                                               C3, var=c(var1, var2), 
                                               studyType = studyType, ...)),
              simplify = FALSE, USE.NAMES = TRUE) 
  
  return(Tt)
}

#------------- RX_VolcanoPlot
RX_VolcanoPlot <- function(DF,
                           logFC_lim = 1,
                           padj_lim = 0.05, 
                           FDRtransform = log10,
                           colors = c("Significative logFC > 0" = "#CB326D", 
                                      "Significative logFC < 0" = "#00a0ab",
                                      "Significative" = "gray24",
                                      "No significative" = "grey"),
                           pval_col = "adj.P.Val",
                           logFC_col = "logFC",
                           point_size = 1,
                           point_shape = 8,
                           label = "ENTREZID",
                           legendTitle = "Significance",
                           ylab = "- log10 (FDR)",
                           xlab = "logFC"){
  # ****************************************************** 
  # Function: RX_VolcanoPlot
  #
  #   This function creates a VolcanoPlot that shows the relationship between 
  #   the log of fold change (logFC) and the negative log of the adjusted p-value 
  #   (FDR) for a dataset. Points are colored based on their significance and 
  #   can be labeled with specific identifiers.
  #
  # Input:
  # -----
  #   - DF (data.frame):
  #         The dataset containing the necessary columns for the plot, including 
  #         logFC, adjusted p-value, and possibly identifiers for labeling points.
  #   - logFC_lim (numeric value, default 1):
  #         The threshold used to define significant points based on logFC.
  #   - padj_lim (numeric value, default 0.05):
  #         The threshold used to define significant points based on the adjusted p-value.
  #   - FDRtransform (function, default log10):
  #         The function used to transform the adjusted p-value on the y-axis of the plot.
  #   - colors (vector of colors, default 
  #             c("Significative logFC > 0" = "#CB326D", 
  #               "Significative logFC < 0" = "#00a0ab",
  #               "Significative" = "gray24",
  #               "No significative" = "grey")):
  #         The colors used to represent different groups of points on the plot.
  #   - pval_col (column name, default "adj.P.Val"):
  #         The name of the column containing the adjusted p-values in the dataset.
  #   - logFC_col (column name, default "logFC"):
  #         The name of the column containing the logFC values in the dataset.
  #   - point_size (numeric value, default 1):
  #         The size of the points on the plot.
  #   - point_shape (numeric value, default 8):
  #         The shape of the points on the plot.
  #   - label (column name, default "ENTREZID"):
  #         The name of the column containing the identifiers used to label the points.
  #   - legendTitle (text, default "Significance"):
  #         The title of the plot's legend.
  #   - ylab (text, default "- log10 (FDR)"):
  #         Label for the y-axis of the plot.
  #   - xlab (text, default "logFC"):
  #         Label for the x-axis of the plot.
  #
  # Output:
  # -----
  #   - (ggplot object):
  #       Volcano plot created with ggplot.
  #
  # Dependencies:
  # -----
  #   - library(ggplot2)
  #   - library(ggrepel)
  #
  # ******************************************************
  
  # We assessed significance
  DF$Sig =  factor(sapply(seq(nrow(DF)), 
                          function(i)
                            ifelse(DF[i,pval_col]>padj_lim, "No significative",
                                   ifelse(DF[i,logFC_col] < -logFC_lim, 
                                          "Significative logFC < 0",
                                          ifelse(DF[i,logFC_col] > logFC_lim,
                                                 "Significative logFC > 0",
                                                 "Significative" )))
                          ))
  
  # We prepare the data
  DF_plot <- DF
  DF_plot$logFDR <- -FDRtransform(DF[,pval_col])
  
  VolcanoPlot <- ggplot(data = DF_plot, aes(x = DF[,logFC_col], y = logFDR)) + 
    # Horizontal bar
    geom_vline(xintercept = c(-logFC_lim, logFC_lim),
               col = "gray", linetype = 'dashed') +
    # Vertical bar
    geom_hline(yintercept = -FDRtransform(padj_lim),
               col = "gray", linetype = 'dashed') + 
    # Points
    geom_point(size = 1, aes(col = Sig), shape=point_shape) + 
    # Fix xlim position
    coord_cartesian(#ylim = c(0, 3),  xlim = c(-1, 1)
      xlim = c(-max(abs(DF[,logFC_col])), max(abs(DF[,logFC_col])))) + 
    # Colors
    scale_color_manual(values = colors)
  
  # Add labels to points
  if(! is.null(label)){
    VolcanoPlot <- VolcanoPlot + 
      # Asign names
      geom_text_repel(aes(DF[,logFC_col], logFDR),
                      label = ifelse(DF$Sig %in% c("Significative logFC < 0",
                                                   "Significative logFC > 0"), 
                                     as.character(DF$SYMBOL),"")) }
  # Theme
  VolcanoPlot <- VolcanoPlot + 
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill = "transparent"),
          axis.line = element_line(colour = "black") ) +
    guides(color =guide_legend(title = legendTitle))
  
  # Labels
  VolcanoPlot <- VolcanoPlot + 
    ylab(ylab) +
    xlab(xlab)  #+ theme(legend.position="none")
  
  return(VolcanoPlot)
}


# ~~~~~~~~~~~~ 05MetaAnalysis ~~~~~~~~~~~~ #

#------------- RX_getSE
RX_getSE <- function(Data 
                     ){
  # ******************************************************** 
  # Function: RX_getSE
  #
  #   This function calculates the standard error (SE) of the estimated 
  #   coefficients in a limma fit object using two methods.
  #
  # Input:
  # -----
  #   - Data (MArrayLM object):
  #         The limma fit object containing the estimated coefficients and the 
  #         necessary information to calculate the SE.
  #
  # Output:
  # -----
  #   - (data.frame):
  #       A data frame containing the estimated coefficients and their 
  #       respective SE.
  #
  # Dependencies:
  # -----
  #   - library(limma)
  #
  # ********************************************************
  
  # Mathod 1:
  print(class(Data))
  TopTab <- topTable(Data, number = "all", confint=TRUE,
                     adjust.method = "fdr", sort.by = "none")
  # SE
  TopTab[, "SE"] <- (TopTab[, "CI.R"] - TopTab[, "CI.L"])/ 3.92
  
  # Mathod 2:
  coefSE = sqrt(Data$s2.post) * Data$stdev.unscaled
  TopTab[, c("coef", "coefSE")] = c(Data$coefficients, coefSE)
  
  return(TopTab)
}

#------------- RX_MetaAnalysis_prep
RX_MetaAnalysis_prep <- function(Studies, # Studies accesion to include
                                 Tissue = c("SAT", "VAT", "LSAT"), # One or several tissues to include in the meta-analysis to include
                                 Contrast, # Contrast to use
                                 Datas, # data
                                 SE_coef = c("SE", "coefSE") # SE to take
){
  # ******************************************************** 
  # Function: RX_MetaAnalysis_prep
  #
  #   This function prepares the data for conducting a gene expression meta-analysis
  #   using the results from multiple studies. It extracts the data for the specified
  #   tissues and contrasts, calculates the standard error (SE) corresponding to the 
  #   estimated coefficients, and creates a logFC matrix and an SE matrix for use in 
  #   the meta-analysis.
  #
  # Input:
  # -----
  #   - Studies (character vector):
  #         The names of the studies included in the meta-analysis.
  #   - Tissue (character vector, one or more of c("SAT", "VAT", "LSAT")):
  #         The tissues to include in the meta-analysis. It can include one or multiple tissues.
  #   - Contrast (character vector):
  #         The contrasts to be used in the meta-analysis.
  #   - Datas (list):
  #         Gene expression data for each study and tissue, organized in a list.
  #   - SE_coef (character vector, one of c("SE", "coefSE")):
  #         The method for obtaining the standard error (SE) in each study's data.
  #
  # Output:
  # -----
  #   - (list):
  #       A list containing the logFC matrix, the SE matrix, and annotation information
  #       for the meta-analysis.
  #
  # Dependencies:
  # -----
  #   - library(dplyr)
  #   - library(limma)
  #
  # ********************************************************
  
  # Check if the arguments are suitable
  match.arg(Tissue, several.ok = TRUE)
  #match.arg(Disease, several.ok = FALSE)
  match.arg(SE_coef, several.ok = FALSE)
  
  # Extract the information from the indicated tissues and contrast
  Data0 = sapply(Studies,
                 function(st) sapply(Tissue,
                                     function(ts) sapply(Contrast,
                                                         function(c) Datas[[st]][[ts]][[c]], 
                                                         simplify = FALSE),
                                     simplify = FALSE),
                 simplify = FALSE)
  
  Data = unlist(unlist(Data0, recursive = FALSE), recursive = FALSE)
  names(Data) = Studies
  Data = Filter(Negate(is.null), Data)
  
  #------------- 1. Preparing input for meta-analysis
  # Get SE
  DataSE = sapply(names(Data), function(x) RX_getSE(Data[[x]]),
                  simplify = FALSE) 
  # Get Gene and annotation
  IDs = unlist(sapply(DataSE, function(x) rownames(x),
                      simplify = TRUE),
               use.names = FALSE)
  SYMBOL = sapply(Data, function(x) 
    x$TopTab[,c("ENTREZID", "ENSEMBL", "SYMBOL")], 
    simplify = FALSE)
  ANOT = do.call("rbind", SYMBOL)
  cat(glue("\n{length(unique(IDs))} of the {length(IDs)} genes have been taken."), "\n") 
  IDs = unique(IDs)
  ANOT = ANOT[match(unique(IDs), ANOT$ENTREZID),]
  rownames(ANOT) = ANOT$ENTREZID
  
  # LogFC matrix
  MatLogFC <- sapply(DataSE, function(x) x[IDs, "logFC"]) 
  rownames(MatLogFC) <- IDs
  
  # SD matrix
  MatSE <- sapply(DataSE, function(x) x[IDs, SE_coef]) 
  rownames(MatSE) <- IDs
  
  # Filter keeping the genes that appear in more than one study
  keep = which(rowSums(!is.na(MatSE)) > 1)
  MatLogFC = MatLogFC[keep,]
  MatSE = MatSE[keep,] 
  
  return(list("MatLogFC" = MatLogFC, "MatSE" = MatSE, Anot = ANOT))
}

