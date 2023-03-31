#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# Functions.R
#
#
# Collection of udeid in this Metha-analysis
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
RX_coind <- function(dict, # Dictionary
                     lst, # List to change based on dictionary information
                     names = FALSE # Keep original names, by default is false
                     ){
  # Map matches between a dictionary and a list
  lst = as.vector(lst)
  out = (sapply(lst, function(x) dict[[x]], USE.NAMES = names))
  return(out)
}

#------------- RX_ifelse
RX_ifelse <- function(cond, # conditional instance
                      true, # Function or object to return if the condition is true
                      false # Function or object to return if the condition is false
                      ){
  if (cond){
    out = true
  }else{
    out = false
  }
  return(out)
}


# ~~~~~~~~~~~~ 01Download ~~~~~~~~~~~~ #

#------------- RX_GetDataGEO
RX_GetDataGEO <- function(StudyAcc, # GEO Accession
                          dir=getwd() # Working directory, by default the current directory is taken
                          ){
  "library(GEOquery)"
  setwd(dir) 
  print(glue("Processing {StudyAcc}:"))
  # Check if the file exists
  if(
    class(
      try(
        {getGEOSuppFiles(StudyAcc, makeDirectory =TRUE,
                         baseDir = dir, fetch_files = TRUE)
          # Set working directory
          setwd(glue("{dir}/{StudyAcc}"))}
      )) == "try-error"){
    # Error
    print(glue("{StudyAcc} could not be downloaded."))
  }else{
    # Decompress files
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
    # Delete the compressed files
    system(paste("rm", files, sep = ' '))
    
    # Download study information
    dir.create("01RawData")
    metadata = getGEO(GEO = StudyAcc)[[1]]
    save(metadata, file = glue("{StudyAcc}_metadata.RData"))
    system( "mv `ls | grep -v 01RawData` 01RawData/")
    
    # Completed
    print(glue("Study {StudyAcc} processed successfully.."))
  }
}

#------------- RX_GetDataArrEx
RX_GetDataArrEx <- function(StudyAcc, # ArrayExpress Accession
                            dir=getwd() # Working directory, by default the current directory is taken
                            ){
  # library(ArrayExpress)
  setwd(dir)
  # Creating the directory to store the data.
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
    # Save the ExpressionSet
    save(exp_set, file = glue("{path0}/{StudyAcc}_raw.rda"))
  }
}


# ~~~~~~~~~~~~ 02Construct ~~~~~~~~~~~~ #

#------------- RX_probe2gene
RX_probe2gene <- function(Data, # ExpressionSet or SummarizedExperiment from microarray study
                          anotdb = org.Hs.eg.db, # Pre-Loaded Annotation Package
                          from = "PROBEID", # Current identifiers
                          to = "ENTREZID", # Identifiers to take
                          multi.to = "first", # Treatment of multiple matches with destination identifiers
                          multi.from = median # Treatment of multiple matches with source identifiers. It corresponds to a specific function, by default the median is taken.
                          ){  
  # library(AnnotationDbi)
  if (class(Data) == "SummarizedExperiment"){ 
    DataOut = Data
    exData = assay(Data)
  }else{
    DataOut = Data 
    exData = exprs(Data)
    #rowData(DataOut) = data.frame(rowData(DataOut))
  }
  # Annotation
  Anot0 = AnnotationDbi::select(x = anotdb,
                                keys = rownames(DataOut),
                                columns = to,
                                keytype = from)
  
  
  # If there is more than one gene per probe, we take the first one.
  ind1 = match(rownames(DataOut), Anot0[,from])
  Anot = Anot0[ind1,]
  
  # Delete probes that do not correspond to any gene. 
  Anot = Anot[which(! is.na(Anot[,to])),]
  
  # Join the annotation with the input data
  Exp0 = data.frame(from = rownames(exData), exData)
  Exp1 = merge(Anot, Exp0, by.x = glue("{from}"), by.y = 1)
  
  # Take the samples name
  Samp = colnames(DataOut)
  
  # If we have more than one probe per gene we use the 'multi.from' function.
  Exp2 = aggregate(Exp1[which(colnames(Exp1) %in% Samp)], 
                   by = list(X= Exp1[,to]),
                   multi.from)
  
  # We have in X the genes according to the indicated nomenclature.
  # Order as in the original,
  ord = match(unique(Exp1[,to]),  Exp2$X)
  Exp3 <- Exp2[ord,]
  ind2 = match(Exp3[,"X"], Anot[,to])
  Anot2 = Anot[ind2,] # Anot 2 es la anotacion final
  
  # Match the data
  ind3 = match(rownames(DataOut), Anot2[,from])
  # Assign the values to the original 
  if (class(Data) == "SummarizedExperiment"){ 
    rowData(DataOut) = Anot2[ind3,]
    DataOut <- DataOut[which(!is.na(rowData(DataOut)$PROBEID),)]
    rownames(DataOut) = Exp3$X
    rownames(Exp3) = Exp3$X
    # We assign the values to the array
    assay(DataOut, withDimnames=FALSE) <- as.matrix(Exp3[,-1])
  }else{
    fData(DataOut) = (Anot2[ind3,])
    DataOut <- DataOut[which(!is.na(fData(DataOut)$PROBEID),)]
    rownames(DataOut) = Exp3$X
    rownames(Exp3) = Exp3$X
    # We assign the values to the array
    exprs(DataOut) <- as.matrix(Exp3[,-1]) 
  }
  # Return
  return(DataOut)
}

#------------- RX_annot
RX_annot <- function(Data, # ExpressionSet or SummarizedExperiment 
                     db = org.Hs.eg.db,  # Pre-Loaded Annotation Package
                     fromAnot="ENTREZID",  # Current identifiers
                     toAnot = c("ENSEMBL","SYMBOL") # Identifiers to take
){
  #library("AnnotationDbi")
  #library("BiocGenerics")
  #library(org.Hs.eg.db)
  Dataout = Data
  # Search for matches
  anot = AnnotationDbi::select(x = db,
                               keys = rownames(Dataout),
                               columns = toAnot,
                               keytype = fromAnot,
                               multiVals = "first")
  
  # Check and deal with multiple matches
  if(! isTRUE(all.equal(rownames(Dataout),anot[,fromAnot]))){  
    ind = BiocGenerics::match(rownames(Dataout), anot[,fromAnot]) 
    if (class(Data) == "SummarizedExperiment"){
      # Combine with the above information
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
# Barplot
RX_resumeBarPlot <- function(Data, # ExpressionSet or SummarizedExperiment 
                             Fac_1, # Name of the first factor (groups the samples on the x-axis and defines the color).
                             Fac_2 = NULL, # Name of the second factor (chance the alpha).
                             Fac_wrap = NULL, # Name of the third factor (group de data).
                             OneBar = TRUE, # Determines how the samples are grouped. IF it is TRUE all the levels of factor 2 are shown in the same bar, if it is FALSE there is one bar per factor.
                             expandY = NULL, # Expansion factor the y-axis (Ex: 0.2 is 20% more).
                             colors = NULL, # Color dictionary for factor 1
                             alpha = c(0.2, 0.6),# Alpha values that make up factor 2
                             ylab = "Número de muestras", # Y label
                             xlab = Fac_1, # X label
                             Fac_1_tit = Fac_1, # Factor one legend title, by default is the name of the variable
                             Fac_2_tit = Fac_2, # Factor two legend title, by default is the name of the variable
                             sep_width=0.4, # Bar spacing
                             numLabels = TRUE, # Add label with the number of samples per group
                             porcentLabels = TRUE, # Add label with the percentage of samples per group with respect to factor 1
                             numpos = 1.6, # Position of the number labels on the y-axis with respect to the end of the bar.
                             porcentpos = 4 # Position of the percentage labels on the y-axis with respect to the end of the bar.
){
  #library(plyr)
  #library(dplyr)
  #library(ggplot)
  # Convertimos los datos
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
    '    BarPlot <- ggplot(Data, aes(x = RX_ifelse(is.null(Fac_1), "",Data[,Fac_1]),  
                                fill= RX_ifelse(is.null(Fac_1), "",Data[,Fac_1]), 
                                color= RX_ifelse(is.null(Fac_1), "",Data[,Fac_1]),
                                alpha= RX_ifelse(is.null(Fac_2), NULL, Data[,Fac_2]))) + 
      geom_bar() '
    # Plot
    BarPlot <- ggplot(Tabla, aes(x = RX_ifelse(is.null(Fac_1), "",Fac_1),
                                 y=Total,
                                 fill= RX_ifelse(is.null(Fac_1), "",Fac_1),
                                 color= RX_ifelse(is.null(Fac_1), "",Fac_1),
                                 alpha= RX_ifelse(is.null(Fac_2), NULL, Fac_2))) + 
      geom_bar(width = sep_width,stat="identity", 
               position = position_stack(reverse = TRUE),
               size  = 1.5
               #color = "#525152" 
               #position = position_dodge()
      ) 
    # Agrupamos
    '    if(! is.null(Fac_wrap)){ 
      BarPlot <- BarPlot + facet_wrap(~Data[, Fac_wrap])} 
    
      # Redimensionamos el eje y a un 10 por ciento más
    if(! is.null(Fac_2) & !is.null(expandY)){ 
      ymax0 = max(table(Data[,Fac_2]))
      ymax = ymax0 +  ymax0 * expandY
      # aumentamos el eje Y
      BarPlot <- BarPlot + ylim(c(0,ymax)) 
    }'
    
    # Etiquetas
    if(isTRUE(numLabels)){
      BarPlot <- BarPlot + 
        # Total
        geom_text(data = Tabla, aes( y = Ypos, label=glue("n = {Total}"), ), 
                  alpha = 1,
                  vjust= numpos,
                  color="black",
                  size=3.0) 
      if(isTRUE(porcentLabels)){
        BarPlot <- BarPlot +
          # Porcentaje
          geom_text(data = Tabla, aes( y = Ypos, label=glue("({Porcentaje} %)"), ), 
                    alpha = 1,
                    vjust=porcentpos,
                    color="black",
                    size=2.5) } 
    }
    # Redimensionamos el eje y a un 10 por ciento más
    if(!is.null(expandY)){ 
      ymax0 = max(Tabla$Ypos)
      ymax = ymax0 +  ymax0 * expandY
      # aumentamos el eje Y
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
               size  = 1.5
      )     
    # Etiquetas
    if(isTRUE(Labels)){
      BarPlot <- BarPlot + 
        # Total
        geom_text(aes(label= glue("n = {Total}"), y = Total), alpha = 1,vjust=1.6,                        
                  color="black", 
                  #hjust=0.5,
                  # define text position and size
                  position = position_dodge2(sep_width),  
                  angle=0, 
                  size=3.0) + 
        # Porcenteje
        geom_text(aes(label= glue("({Porcentaje} %)"), y = Total), alpha = 1,vjust=4,                        
                  color="black", 
                  #hjust=0.5,
                  # define text position and size
                  position = position_dodge2(sep_width),  
                  angle=0, 
                  size=2.5) }
    # Redimensionamos el eje y a un 10 por ciento más
    if(!is.null(expandY)){ 
      ymax0 = max(Tabla$Total)
      ymax = ymax0 +  ymax0 * expandY
      # aumentamos el eje Y
      BarPlot <- BarPlot + ylim(c(0,ymax)) 
    }   
  }  
  # Agrupamos
  if(! is.null(Fac_wrap)){ 
    BarPlot <- BarPlot + facet_wrap(~Fac_wrap)}
  
  # Ajuste   
  # COlores 
  if(! is.null(colors)){
    BarPlot <- BarPlot + 
      # Asignamos colores
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
    # Nombres ejes
    ylab(ylab) +
    xlab(xlab) + 
    # Leyenda 
    guides(fill = guide_legend(title=Fac_1_tit,
                               order = 1),
           alpha = guide_legend(title=Fac_2_tit,
                                order = 2),
           color = "none")
  
  return(BarPlot)
}
