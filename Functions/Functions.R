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
  # Función: RX_coind()       
  #     
  #   Esta función toma una lista y la modifica según la información de un diccionario.
  #
  #   Input:
  #   -----
  #    - dict (vector con nombres):
  #         Diccionario donde las claves corresponden a los elementos en 
  #         la lista `lst`, y los valores son las modificaciones
  #         deseadas para esos elementos.
  #    - lst (lista):
  #         Lista que se modificará según la información del diccionario.
  #    - names (lógico, opcional, por defecto falso):  
  #         Determina si se mantendrán los nombres originales de los elementos
  #         de la lista. .
  # 
  #   Output:
  #   -----
  #   - (vector) Lista modificada según la información del diccionario.
  # 
  # Dependencias:
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
  # Función: RX_ifelse()       
  #     
  #   Esta función implementa una versión simplificada de la función ifelse en R.
  #
  #   Input:
  #   -----
  #    - cond: 
  #         Condición lógica que determina si se debe devolver el resultado
  #         verdadero o falso.
  #    - true: 
  #         Valor a devolver si la condición es verdadera.
  #    - false:
  #         Valor a devolver si la condición es falsa.
  # 
  #   Output:
  #   -----
  #   - (lógico) Valor correspondiente a la condición evaluada.
  # 
  # Dependencias:
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
  # Función: RX_GetDataGEO
  #     
  #   Esta función se utiliza para obtener datos del Gene Expression Omnibus
  #   (GEO) mediante el acceso GEO proporcionado. Descarga los archivos
  #   necesarios, los descomprime y guarda la información del estudio en un
  #   archivo de metadatos.
  #
  # Input:
  # -----
  #   - StudyAcc (character): 
  #         El identificador de acceso GEO del estudio.
  #   - dir (character): 
  #         Directorio de trabajo. Por defecto, se toma el directorio actual.
  # Output:
  # -----
  #   none
  #
  # Dependencias:
  # -----
  #   library(GEOquery)
  #
  # ******************************************************
  
  setwd(dir) 
  print(glue("Procesando {StudyAcc}:"))
  
  # Verificar si el archivo existe
  if(
    class(
      try(
        {
          getGEOSuppFiles(StudyAcc, makeDirectory = TRUE,
                          baseDir = dir, fetch_files = TRUE)
          # Establecer el directorio de trabajo
          setwd(glue("{dir}/{StudyAcc}"))
        }
      )) == "try-error"
  ){
    # Error
    print(glue("{StudyAcc} no pudo ser descargado."))
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
    # Eliminar los archivos comprimidos
    system(paste("rm", files, sep = ' '))
    
    # Descargar información del estudio
    dir.create("01RawData")
    metadata = getGEO(GEO = StudyAcc)[[1]]
    save(metadata, file = glue("{StudyAcc}_metadata.RData"))
    system( "mv `ls | grep -v 01RawData` 01RawData/")
    
    # Completado
    print(glue("Estudio {StudyAcc} procesado exitosamente."))
  }
}

  
#------------- RX_GetDataArrEx
RX_GetDataArrEx <- function(StudyAcc, dir = getwd()){
  # ****************************************************** 
  # Función: RX_GetDataArrEx
  #     
  #   Esta función se utiliza para obtener datos del ArrayExpress mediante el
  #   acceso ArrayExpress proporcionado. Crea un directorio para almacenar los
  #   datos, descarga el conjunto de datos y guarda el objeto ExpressionSet en
  #   un archivo .rda.
  #
  # Input:
  # -----
  #   - StudyAcc (character): 
  #         El identificador de acceso ArrayExpress del estudio.
  #   - dir (character): 
  #         Directorio de trabajo. Por defecto, se toma el directorio actual.
  # Output:
  # -----
  #   none
  #
  # Dependencias:
  # -----
  #   library(ArrayExpress)
  #
  # ****************************************************** 
  
  setwd(dir)
  # Crear el directorio para almacenar los datos
  path0 = glue("{StudyAcc}/01RawData")
  system(glue("mkdir {StudyAcc} {path0}"))
  
  # Descarga de datos
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
  # Función: RX_probe2gene
  #     
  #   Esta función se utiliza para anotar los nombres de genes correspondientes
  #   a las sondas o "probes" en un objeto ExpressionSet o SummarizedExperiment
  #   de un estudio de microarrays. para realizar la asignación utiliza una 
  #   base de datos de anotaciones previmanete cargada.
  #
  # Input:
  # -----
  #   - Data (ExpressionSet o SummarizedExperiment):
  #         Objeto ExpressionSet o SummarizedExperiment del estudio.
  #   - anot_tb (data.frame, opcional): 
  #         Tabla de anotaciones. Si se proporciona, se utiliza en lugar
  #         de la base de datos de anotaciones pre-cargada.
  #   - anotdb (AnnotationDbi, por defecto "org.Hs.eg.db".): 
  #         Paquete de anotaciones pre-cargado. 
  #   - from (character, or defecto "PROBEID"): 
  #         Identificadores actuales de las sondas/probes. 
  #   - to (character, por defecto "ENTREZID".): 
  #         Identificadores a los que se asignarán los nombres de genes.
  #         
  #   - multi.to (character, por defecto, se toma el primer valor encontrado): 
  #         Tratamiento de las coincidencias múltiples con los identificadores
  #         de destino. 
  #   - multi.from (function, por defecto se utiliza la función "median"): 
  #         Tratamiento de las coincidencias múltiples con los identificadores
  #         de origen. .
  # Output:
  # -----
  #   - ExpressionSet o SummarizedExperiment: 
  #       Objeto modificado con los nombres de genes asignados a las sondas.
  #
  # Dependencias:
  # -----
  #   library(AnnotationDbi)
  #
  # ****************************************************** 
  
  # Comprobación del tipo de objeto
  if (class(Data) == "SummarizedExperiment"){ 
    DataOut = Data
    exData = assay(Data)
  }else{
    DataOut = Data 
    exData = exprs(Data)
    #rowData(DataOut) = data.frame(rowData(DataOut))
  }
  # Anotacion
  if(is.null(anot_tb)){
    Anot0 = AnnotationDbi::select(x = anotdb,
                                  keys = rownames(DataOut),
                                  columns = to,
                                  keytype = from)
  }else{
    Anot0 = anot_tb
  }
  
  # Si hay más de un gen por sonda, se toma el primero.
  ind1 = match(rownames(DataOut), Anot0[,from])
  Anot = Anot0[ind1,]
  
  # Eliminar las sondas que no corresponden a ningún gen. 
  Anot = Anot[which(! is.na(Anot[,to])),]
  
  # Unir la anotación con los datos de entrada
  Exp0 = data.frame(from = rownames(exData), exData)
  Exp1 = merge(Anot, Exp0, by.x = glue("{from}"), by.y = 1)
  
  # Tomar el nombre de las muestras
  Samp = colnames(DataOut)
  
  # Si hay más de una sonda por gen usar la función 'multi.from'.
  Exp2 = aggregate(Exp1[which(colnames(Exp1) %in% Samp)], 
                   by = list(X= Exp1[,to]),
                   multi.from)
  
  # Tenemos en X los genes según la nomenclatura indicada.
  # Ordenar como en el original.
  ord = match(unique(Exp1[,to]),  Exp2$X)
  Exp3 <- Exp2[ord,]
  ind2 = match(Exp3[,"X"], Anot[,to])
  Anot2 = Anot[ind2,] # Anot 2 es la anotacion final
  
  # Asignar los valores a la matriz original
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
  # Devolver el objeto modificado
  return(DataOut)
}

#------------- RX_annot
RX_annot <- function(Data,  
                     db = org.Hs.eg.db,  
                     fromAnot="ENTREZID", 
                     toAnot = c("ENSEMBL","SYMBOL") 
                     ){
  # ****************************************************** 
  # Función: RX_annot
  #     
  #   Esta función se utiliza para realizar la anotación de identificadores en 
  #   un objeto de tipo ExpressionSet o SummarizedExperiment. Busca 
  #   coincidencias entre los identificadores actuales y los identificadores 
  #   deseados utilizando una base de datos de anotación pre-cargada.
  #
  # Input:
  # -----
  #   - Data (ExpressionSet o SummarizedExperiment): 
  #         Objeto que contiene los datos a anotar.
  #   - db (AnnotationDbPackage): 
  #         Paquete de base de datos de anotación pre-cargado.
  #   - fromAnot (character, por defecto, es "ENTREZID"): 
  #         Identificador actual de los datos.
  #   - toAnot (character vector, por defecto, es "ENSEMBL" y "SYMBOL"): 
  #         Identificadores deseados a tomar. 
  #
  # Output:
  # -----
  #   - (ExpressionSet o SummarizedExperiment):
  #         Objeto de datos anotado con los identificadores deseados.
  #
  # Dependencias:
  # -----
  #   library(AnnotationDbi)
  #   library(BiocGenerics)
  #   library(org.Hs.eg.db)
  #
  # ****************************************************** 
  
  Dataout = Data
  # Buscar coincidencias
  anot = AnnotationDbi::select(x = db,
                               keys = rownames(Dataout),
                               columns = toAnot,
                               keytype = fromAnot,
                               multiVals = "first")
  
  # Verificar y manejar coincidencias múltiples
  if(! isTRUE(all.equal(rownames(Dataout),anot[,fromAnot]))){  
    ind = BiocGenerics::match(rownames(Dataout), anot[,fromAnot]) 
    if (class(Data) == "SummarizedExperiment"){
      # Combinar con la información anterior
      rowData(Dataout) =cbind(anot[ind,],
                              rowData(Dataout)[which(!colnames(rowData(Dataout))
                              %in% c(fromAnot,toAnot))])
    }else{
      # Combinar con la información anterior
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
                             ylab = "Número de muestras",
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
  # Función: RX_resumeBarPlot
  #     
  #   Esta función se utiliza para generar un gráfico de barras que resume la
  #   información en un objeto de tipo ExpressionSet o SummarizedExperiment.
  #   El gráfico muestra la distribución de muestras en grupos definidos por
  #   dos factores (Fac_1 y Fac_2) y, opcionalmente, un tercer factor (Fac_wrap).
  #
  # Input:
  # -----
  #   - Data (ExpressionSet o SummarizedExperiment): 
  #         Objeto que contiene los datos a visualizar.
  #   - Fac_1 (character): 
  #         Nombre del primer factor que agrupa las muestras en el eje x,
  #         y define el color de las barras.
  #   - Fac_2 (character, opcional): 
  #         Nombre del segundo factor que cambia la transparencia de las barras.
  #   - Fac_wrap (character, opcional): 
  #         Nombre del tercer factor que agrupa los datos y crea subgráficos.
  #   - OneBar (logical, por defecto TRUE): 
  #         Determina cómo se agrupan las muestras.
  #         Si es TRUE, todos los niveles del factor 2 se muestran en una misma
  #         barra; si es FALSE, hay una barra por cada nivel del factor.
  #   - expandY (numeric, opcional): 
  #         Factor de expansión del eje y (por ejemplo, 0.2
  #         representa un 20% más).
  #   - colors (vector, opcional): 
  #         Diccionario de colores para el factor 1.
  #   - alpha (numeric vector, por defecto c(0.2, 0.6)): 
  #         Valores alpha que conforman el factor 2 
  #         (transparencia de las barras).
  #   - ylab (character, por defecto "Número de muestras"): 
  #         Etiqueta del eje y.
  #   - xlab (character, por defecto Fac_1): 
  #         Etiqueta del eje x.
  #   - Fac_1_tit (character, por defecto Fac_1): 
  #         Título de la leyenda del factor 1.
  #   - Fac_2_tit (character, por defecto Fac_2): 
  #         Título de la leyenda del factor 2.
  #   - sep_width (numeric, por defecto 0.4): 
  #         Espaciado entre barras.
  #   - Labels (logical, por defecto TRUE): 
  #         Determina si se agregan etiquetas en las barras.
  #   - numLabels (logical, por defecto TRUE): 
  #         Agrega etiquetas con el número de muestras por grupo.
  #   - porcentLabels (logical, por defecto TRUE): 
  #         Agrega etiquetas con el porcentaje de muestras por grupo en
  #         relación con el factor 1.
  #   - numpos (numeric, por defecto 1.6): 
  #         Posición de las etiquetas numéricas en el eje y con respecto al
  #         final de la barra.
  #   - porcentpos (numeric, por defecto 4): 
  #         Posición de las etiquetas de porcentaje en el eje y con respecto al
  #         final de la barra.
  #   - size (numeric, por defecto 1): 
  #         Tamaño del borde de las barras.
  #   - text_size (numeric, por defecto 3.5): 
  #         Tamaño del texto principal en el gráfico.
  #   - porcent_size (numeric, por defecto 3.5): 
  #         Tamaño del porcentaje sobre el gráfico.
  # Output:
  # -----
  #   - (objeto de la clase "ggplot"):
  #         Gráfico de barras obtenido con ggplot.
  #
  # Dependencias:
  # -----
  #   library(plyr)
  #   library(dplyr)
  #   library(ggplot2)
  #
  # ****************************************************** 
  
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

  Tabla <- Tabla %>%
    ddply(.variables = c("Fac_1", RX_ifelse(is.null(Fac_wrap), NULL,
                                            "Fac_wrap")), 
          transform, Porcentaje = round(Total / sum(Total) * 100, 1))
  
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
    
    # Etiquetas
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
          # Porcentaje
          geom_text(data = Tabla, aes( y = Ypos,
                                       label=glue("({Porcentaje} %)"), ), 
                    alpha = 1,
                    vjust=porcentpos,
                    color="black",
                    size=porcent_size) } 
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
               size  = size
      )     
    # Etiquetas
    if(isTRUE(Labels)){
      BarPlot <- BarPlot + 
        # Total
        geom_text(aes(label= glue("N = {Total}"), y = Total),
                  alpha = 1,vjust=1.6,                        
                  color="black", 
                  #hjust=0.5,
                  # define text position and size
                  position = position_dodge2(sep_width),  
                  angle=0, 
                  size=text_size) + 
        # Porcenteje
        geom_text(aes(label= glue("({Porcentaje} %)"), y = Total),
                  alpha = 1,vjust=4,                        
                  color="black", 
                  #hjust=0.5,
                  # Define el tamaño y posición del texto
                  position = position_dodge2(sep_width),  
                  angle=0, 
                  size=porcent_size) }
    # Redimensionar el eje y a un 10 por ciento más
    if(!is.null(expandY)){ 
      ymax0 = max(Tabla$Total)
      ymax = ymax0 +  ymax0 * expandY
      # Aumentar el eje Y
      BarPlot <- BarPlot + ylim(c(0,ymax)) 
    }   
  }  
  # Agrupar
  if(! is.null(Fac_wrap)){ 
    BarPlot <- BarPlot + facet_wrap(~Fac_wrap)}
  
  # Ajuste   
  # Colores 
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
    # Nombre de los ejes
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

#------------- RX_dataforPCA 
RX_dataforPCA <- function(data){
  # ****************************************************** 
  # Función: RX_dataforPCA
  #
  #   Esta función realiza el cálculo de PCA sobre los datos proporcionados.
  #
  # Input:
  # -----
  #   - data (matrix o data frame):
  #         Datos sobre los que computar la PCA.
  #
  # Output:
  # -----
  #   - (objeto de clase "prcomp"):
  #         Objeto que contiene los resultados del análisis de PCA.
  #
  # Dependencias:
  # -----
  #   library(stats)
  #
  # ****************************************************** 
  
  # Calculo de la matriz transpuesta
  TData = t(data)
  # Cálculo del PCA 
  DataPCA = prcomp(TData,scale.=TRUE)
  return(DataPCA)
}

#------------- RX_dataforPCA
RX_elimrows0 <- function(data){ 
  # ****************************************************** 
  # Función: RX_elimrows0
  #
  #   Esta función elimina las filas de una matriz o data frame que contengan 
  #   únicamente valores iguales a 0.
  #
  # Input:
  # -----
  #   - data (matrix o data frame):
  #         Datos de los cuales a procesar.
  #
  # Output:
  # -----
  #   - (matrix o data frame):
  #         Datos modificados (eliminando filas 0).
  #
  # Dependencias:
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
  # Función: RX_plotPCA
  #
  #   Esta función crea un gráfico de PCA a partir de un objeto PCA y datos 
  #   fenotípicos.
  #
  # Input:
  # -----
  #   - PCAData (objeto PCA):
  #         Objeto PCA obtenido a partir de la función prcomp.
  #   - phenoData (data frame):
  #         Datos fenotípicos utilizados para agregar información a los
  #         puntos del gráfico.
  #   - shade (numérico, por defecto  0.3):
  #         Valor de sombreado para los puntos de fondo. 
  #   - shapes (vector numérico, por defecto c(16,17)):
  #         Vector de códigos de formas para los puntos del gráfico. 
  #   - labels (character, opcional):
  #         Nombre de la variable que se desea poner en las etiquetas de los
  #         puntos. 
  #   - colors (vector de colores, opcional):
  #         Vector de colores para los puntos del gráfico. 
  #   - Fac_wrap (character, opcional):
  #         Nombre de la variable que se utiliza para separar los puntos en 
  #         múltiples paneles. 
  #   - colorColumn (character, opcional):
  #         Nombre de la columna que se utilizará para asignar colores a los
  #         puntos. 
  #   - shapeColumn (character, opcional):
  #         Nombre de la columna que se utilizará para asignar formas a los
  #         puntos. 
  #   - colorLegend (character, por defecto "colorColumn"):
  #         Nombre de la columna que se utilizará para el título de la leyenda
  #         de colores. Si es NULL, no se mostrará la leyenda de colores. 
  #   - shapeLegend (character, por defecto "colorColumn"):
  #         Nombre de la columna que se utilizará para el título de la leyenda
  #         de formas. Si es NULL, no se mostrará la leyenda de formas. 
  #   - alpha (numérico, por defecto 0.8):
  #         Valor de transparencia de los puntos del gráfico. 
  #   - size (numérico, por defecto 5):
  #         Tamaño de los puntos del gráfico.
  #
  # Output:
  # -----
  #   - (objeto de la clase "ggplot"):
  #         Gráfico de PCA obtenido con ggplot.
  #
  # Dependencias:
  # -----
  #   - library(ggplot2)
  #   - library(ggrepel)
  #   - library(stats)
  #
  # ****************************************************** 
  
  # Conversión de los datos 
  inputData <- data.frame(PCAData$x[,c(1,2)], phenoData) 
  
  #  Creación del objeto ggplot
  PCAplot = ggplot(inputData) 
  
  # Si separamos por un factor
  if (! is.null(Fac_wrap)){
    PCAplot = PCAplot +
      # Puntos de fondo
      geom_point(data= inputData[,which(colnames(inputData)!=Fac_wrap)],
                 aes(PC1, PC2,
                     colour  = "grey",
                     shape = Gender,
                     alpha = shade,
                     size = 5))
    } 
  # Agregar puntos del PCA      
  PCAplot = PCAplot +
    geom_point(aes(PC1, PC2,
                     colour  = RX_ifelse(is.null(colorColumn),
                                         NULL,inputData[,colorColumn]),
                     shape = RX_ifelse(is.null(shapeColumn),
                                       NULL,inputData[,shapeColumn])
                   ), 
               alpha = alpha,
               size = size)  + 
      # Tema
      theme(plot.background = element_rect(fill="transparent"),
            panel.background = element_rect(fill = "transparent"),
            axis.line = element_line(colour = "black")) +
      # Colores
      scale_color_manual(values = colors) +
      # Forma
      scale_shape_manual(values = shapes) 
    
    if(! is.null(Fac_wrap)){
      # Division por tejidos
      PCAplot = PCAplot +
        facet_wrap(facets = c(Fac_wrap), ncol = 2, scales = "free")} 
  
      # Titulo leyendas
      # Leyenda (intensidad no se crea sola)
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
      xlab(glue("PCA1 ({round(summary(PCAData)$importance[2,1], 4) *100}% de varianza explicada)"))+
      ylab(glue("PCA2 ({round(summary(PCAData)$importance[2,2], 4) *100}% de varianza explicada)"))
  
    
  # Calculamos la escala
  escalx = max(inputData$PC1) - min(inputData$PC1)
  escaly = max(inputData$PC2) - min(inputData$PC2) 
  # Etiquetas de los puntos
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
  # Función: RX_multiplotPCA
  #
  #   Esta función elimina las filas de una matriz o data frame que contengan 
  #   únicamente valores iguales a 0.
  #
  # Input:
  # -----
  #   - Data (ExpressionSet o SummarizedExperiment):
  #         Datos de expresión para realizar el análisis de PCA.
  #   - multiplot (lógico, por defecto FALSE):
  #         Indica si se debe crear un PCA para cada nivel del factor o
  #         un único PCA para todos los datos.
  #   - factor (character, opcional):
  #         Nombre de la columna que se utilizará para separar los datos
  #         en múltiples PCAs. Solo se usa si multiplot = TRUE.
  #   - ...:
  #         Otros argumentos que se pasan a la función RX_plotPCA.
  #
  # Output:
  # -----
  #   - (objeto de la clase "ggplot" o lista de objetos "ggplot"):
  #         Gráfico(s) de PCA obtenido(s) con ggplot.
  #
  # Dependencias:
  # -----
  #   - library(glue)
  #   - library(ggplot2)
  #   - library(ggrepel)
  #   - library(stats)
  #
  # ****************************************************** 
  if(isFALSE(multiplot)){
    # Creamos un único PCA 
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
    # Creamos un PCA para cada nivel del factor
    # Niveles a tomar
    Fac = glue("Data${factor}")
    Levs = eval(parse(text=Fac))
    # División
    T1 = "Data[,which(Levs == x)]" 
    
    # Martiz
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
    # Datos fenotípicos
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
  # Función: RX_dataforboxplot
  #
  #   Esta función elimina las filas de una matriz o data frame que contengan 
  #   únicamente valores iguales a 0.
  #
  # Input:
  # -----
  #   - data (matriz):
  #         Matriz de expresión.
  #   - transformacion (character, opcional):
  #         Tipo de transformación a aplicar a los datos.
  #         Puede ser "log2" o "log10". Por defecto, es NULL, lo que significa
  #         que no se aplicará ninguna transformación.
  #   - phenoData (data frame o matriz, opcional):
  #         Datos fenotípicos utilizados para agregar información adicional a 
  #         los boxplots.
  #   - ncol (entero, por defecto 1):
  #         Índice de la columna en phenoData que contiene los nombres de
  #         las muestras.
  #
  # Output:
  # -----
  #   - (objeto de la clase "data.frame"):
  #         Datos preparados para crear un gráfico de boxplot.
  #
  # Dependencias:
  # -----
  #   - library(dplyr)
  #   - library(tidyr)
  #
  # ****************************************************** 

  # Aplicar transformación si se especifica
  if(! is.null(transformacion) ){
    if(transformacion == "log2"){
      data = log2(data)
    }else{
      if(transformacion == "log10"){
        data = log10(data)
      }
    }}
  # Convertir  los datos a data.frame
  data = data.frame(data)
  
  # Preparar los datos para la función pivot_longer
  Dt = data %>% pivot_longer(cols = c(colnames(data)),
                             names_to = "Sample",
                             values_to = "exp")
  
  # Ordenar los datos según los originales
  Dt$Sample = factor(Dt$Sample,
                     levels = colnames(data))
  Dt = Dt[(order(Dt$Sample)),]
  
  # Unir los datos fenotípicos si se especifican
  if(! is.null(phenoData)){
    Dt = merge(x = Dt, y= phenoData, 
               by.x = "Sample", by.y = ncol,
               sort = FALSE)}
  
  return(Dt)
}

#------------- RX_BoxPlot 
RX_BoxPlot <- function(Data,
                       Fac_fill = "Group",
                       Fac_alpha = "Gender",
                       Fac_wrap= "Tissue",
                       colors = c("#90C432", "#2494B5"),
                       alpha = c(0.95, 0.45),
                       shape = 19,
                       dir = "h",# 8
                       ncol = 1
){
  # ****************************************************** 
  # Función: RX_BoxPlot
  #
  # Crea un gráfico de boxplot.
  #
  # Input:
  # -----
  #   - Data (data frame):
  #         Datos utilizados para crear el gráfico de boxplot.
  #   - Fac_fill (character, por defecto "Group"):
  #         Nombre de la columna utilizada para el relleno de los boxplots.
  #   - Fac_alpha (character, por defecto "Gender"):
  #         Nombre de la columna utilizada para la transparencia de los
  #         boxplots.
  #   - Fac_wrap (character, por defecto "Tissue"):
  #         Nombre de la columna utilizada para la división de los boxplots
  #         en múltiples paneles.
  #   - colors (vector de colores, por defecto c("#90C432", "#2494B5")):
  #         Colores utilizados para el relleno de los boxplots.
  #   - alpha (vector numérico, por defecto c(0.95, 0.45)):
  #         Valores de transparencia utilizados para los boxplots.
  #   - shape (entero, por defecto 19):
  #         Código de forma utilizado para los puntos atípicos en los boxplots.
  #   - dir (character, por defecto "h"):
  #         Dirección de la división de los paneles. "h" para horizontal y
  #         "v" para vertical.
  #   - ncol (entero, por defecto 1):
  #         Número de columnas utilizado para la disposición de los paneles.
  #
  # Output:
  # -----
  #   - (objeto de la clase "ggplot"):
  #         Boxplot obtenido con ggplot.
  #
  # Dependencias:
  # -----
  #   - library(ggplot2)
  #
  # ****************************************************** 
  
  BoxPlot = (
    # Objeto ggplot
    ggplot(Data,
           aes(x = Sample, y = exp)) +
      # Boxplot
      geom_boxplot(aes(fill = Data[,Fac_fill],
                       alpha = Data[,Fac_alpha]),
                   outlier.colour = "black", outlier.shape = shape,
                   outlier.alpha = 1)+ 
      # Tema
      theme(plot.background = element_rect(fill="transparent"),
            panel.background = element_rect(fill = "transparent"),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 1, size = 8)) +
      # Colores de relleno
      scale_fill_manual(values = colors) +  
      # Intensidad de colors
      scale_alpha_discrete(range = alpha) + 
      # Division por tejidos
      facet_wrap(~get(Fac_wrap), ncol = ncol, scales = "free", dir = dir) +
      # Titulo leyendas
      # Leyenda (intensidad no se crea sola)
      guides(fill = guide_legend(title="Grupo",
                                 order = 1),
             alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),100,0,
                                                           alpha=alpha),
                                                  colour="black"),
                                title = "Sexo",
                                order = 2)) +
      # Titulo ejes
      xlab("Muestras") +
      ylab("log2(expresión)") +
      # Etiquetas 2 niveles eje x
      scale_x_discrete(guide = guide_axis(n.dodge = 2))
  )
  return(BoxPlot)
}


#------------- RX_Clustering 
RX_Clustering <- function(expressData,
                          method = c("euclidean", "correlation")
                          ){
  # ****************************************************** 
  # Función: RX_Clustering
  #
  #   Esta función realiza computa un clustering jerárquico utilizando la
  #   correlación o la distancia euclidiana.
  #
  # Input:
  # -----
  #   - expressData (matriz o data frame):
  #         Datos de expresión utilizados para el clustering.
  #   - method (character, por defecto c("euclidean", "correlation")):
  #         Método utilizado para el cálculo de la distancia.
  #         Puede ser "euclidean" para distancia euclidiana o
  #         "correlation" para el coeficiente de correlación de Pearson
  #
  # Output:
  # -----
  #   - (objeto de la clase "hclust"):
  #         Objeto de clustering jerárquico obtenido.
  #
  # Dependencias:
  # -----
  #   - library(stats)
  #
  # ******************************************************
  
  # Cluster herarquico y trasformación de los datos
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
                           ncol = 1){ # Añade información fenotipica a un dendograma
  # ****************************************************** 
  # Función: RX_datatodendo
  #
  #   Esta función añade información fenotípica a un dendrograma.
  #
  # Input:
  # -----
  #   - hc (objeto de la clase "hclust"):
  #         Objeto de clustering jerárquico obtenido previamente.
  #   - phenoData (data frame, opcional):
  #         Datos fenotípicos utilizados para agregar información a
  #         los nodos del dendrograma.
  #   - ncol (numérico, por defecto 1):
  #         Índice de columna en phenoData que coincide con los nombres 
  #         de los nodos en el objeto hc.
  #
  # Output:
  # -----
  #   - (objeto de la clase "dendro"):
  #         Objeto de dendrograma con la información fenotípica añadida.
  #
  # Dependencias:
  # -----
  #   - library(ggdendro)
  #
  # ******************************************************
  # Convertimos el objeto de clustering jerárquico a dendrograma
  dendoData0 <- as.dendrogram(hc)
  dendoData <- dendro_data(dendoData0)
  
  # Añadimos los datos fenotipicos
  if (!is.null(phenoData)){
    dendoData$labels <- data.frame(merge(dendoData$labels, phenoData, 
                                         by.x = "label", by.y=ncol,
                                         sort = FALSE))
  }
  # Devolvemos el objeto 
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
  # Función: RX_dendoPlot
  #
  #   Esta función dibuja un dendrograma a partir de los datos proporcionados en dendoData.
  #
  # Input:
  # -----
  #   - dendoData (objeto de la clase "dendro"):
  #         Objeto de dendrograma con la información necesaria para dibujar el 
  #         dendrograma.
  #   - dir (carácter, por defecto "lr"):
  #         Dirección de dibujo del dendrograma. Puede ser "lr" (de izquierda 
  #         a derecha), "rl" (de derecha a izquierda), "tb" (de arriba hacia  
  #         abajo)o "bt" (de abajo hacia arriba).
  #   - colorPal (vector de caracteres, opcional):
  #         Paleta de colores utilizada para los puntos en las hojas del 
  #         dendrograma.
  #   - colorOutPal (vector de caracteres, opcional):
  #         Paleta de colores utilizada para las ramas y segmentos del
  #         dendrograma.
  #   - colorOutColumn (carácter, opcional):
  #         Paleta de colores utilizada para las ramas y segmentos del
  #         dendrograma.
  #   - colorOutColumn (carácter, opcional):
  #         Nombre de la columna en dendoData que contiene la variable
  #         utilizada para el color de las ramas y segmentos.
  #   - branchSize (numérico, opcional):
  #         Tamaño de las ramas del dendrograma.
  #   - labelSize (numérico, por defecto 2):
  #         Tamaño de las etiquetas en el dendrograma.
  #   - moveLabel (numérico, por defecto 0.001):
  #         Desplazamiento vertical de las etiquetas en el dendrograma.
  #   - expandY (numérico, por defecto 0.1):
  #         Factor de expansión vertical del dendrograma.
  #   - expandYlbs (numérico, por defecto 0):
  #         Factor de expansión vertical de las etiquetas en relación con
  #         la longitud total del dendrograma.
  #   - colorColumn (carácter, opcional):
  #         Nombre de la columna en dendoData que contiene la variable
  #         utilizada para el color de los puntos en las hojas.
  #   - leaves (lógico, por defecto FALSE):
  #         Indica si se deben mostrar los puntos en las hojas del dendrograma.
  #   - shapeColumn (carácter, opcional):
  #         Nombre de la columna en dendoData que contiene la variable
  #         utilizada para la forma de los puntos en las hojas.
  #   - leavesSize (numérico, por defecto 2):
  #         Tamaño de los puntos en las hojas del dendrograma.
  #   - leavesStroke (numérico, por defecto 1):
  #         Grosor del contorno de los puntos en las hojas del dendrograma.
  #   - leavesShapes (vector de caracteres, opcional):
  #         Formas utilizadas para los puntos en las hojas del dendrograma.
  #   - colorLegend (carácter, opcional):
  #         Nombre de la columna en dendoData que contiene la variable 
  #         utilizada para la leyenda de color.
  #   - colorOutLegend (carácter, opcional):
  #         Nombre de la columna en dendoData que contiene la variable
  #         utilizada para la leyenda de color de las ramas.
  #   - shapeLegend (carácter, opcional):
  #         Nombre de la columna en dendoData que contiene la variable
  #         utilizada para la leyenda de forma.
  #
  # Output:
  # -----
  #   - (objeto de la clase "ggplot"):
  #         Objeto de ggplot que representa el dendrograma.
  #
  # Dependencias:
  # -----
  #   - library(ggnewscale)
  #   - library(ggplot2)
  #
  # ****************************************************** 
  
  # Parametros
  dir <- match.arg(dir)
  ybreaks   <- pretty(segment(dendoData)$y, n = 5)
  ymax      <- max(segment(dendoData)$y)
  
  # Dendrograma básico
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
  
  # Etiquetas
  # Fijamos los parametros
  #{
  nLabs = nrow(dendoData$labels)
  angle <- rep(0, nLabs)
  hjust <- 0
  if (dir %in% c("tb", "bt")) {angle <- angle + 45 }
  if (dir %in% c("tb", "rl")) {hjust <- 1 }
  
  labelParams = list(angle = angle, hjust = hjust, vjust = 0.5)
  dendoData$labels$angle <- labelParams$angle 
  #}
  
  # Añadimos las etiquetas
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
  # Asignamos la paleta de colores
  DendoPlot <- DendoPlot + guides(color = "none")
  
  if (!is.null(colorPal)) {
    DendoPlot <- DendoPlot + scale_color_manual(values = colorOutPal)
  } 
  
  #Añadimos los puntos  
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
  
  # Modificamos la longitud
  ylim <- ymax * expandY#-round(ymax * expandY, 1)
  DendoPlot <- DendoPlot + expand_limits(y = ylim)
  
  # Leyendas
  
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
  # Función: RX_multiplotDendo
  #
  #   Esta función crea dendrogramas a partir de los datos proporcionados y 
  #   los organiza en una o múltiples figuras, dependiendo de los parámetros 
  #   especificados.
  #
  # Input:
  # -----
  #   - Data (objeto de clase "SummarizedExperiment" o "matrix"):
  #         Los datos utilizados para construir los dendrogramas. 
  #   - method (carácter, por defecto "correlation"):
  #         Método utilizado para calcular la matriz de distancias para
  #         el clustering jerárquico. Puede ser "correlation",
  #         "euclidean", "manhattan" u otros métodos disponibles en la
  #         función "dist".
  #   - multiplot (lógico, por defecto FALSE):
  #         Indica si se deben crear dendrogramas individuales para cada 
  #         nivel del factor especificado.
  #   - factor (carácter, opcional):
  #         Nombre del factor utilizado para crear dendrogramas individuales. 
  #         Solo se utiliza si "multiplot" es TRUE.
  #   - wrap (lógico, por defecto FALSE):
  #         Indica si los dendrogramas individuales deben ser envueltos
  #         en un diseño de cuadrícula.
  #   - ... :
  #         Otros argumentos que se pasan a la función "RX_dendoPlot" 
  #         para personalizar los dendrogramas.
  #
  # Output:
  # -----
  #   - (objeto de la clase "ggplot"):
  #         Objeto de ggplot que representa el dendrograma.
  #
  # Dependencias:
  # -----
  #   - library(glue)
  #   - library(ggnewscale)
  #   - library(ggplot2)
  #   - library(stats)
  #
  # ****************************************************** 
  
  if(isFALSE(multiplot)){
    # Creamos un único Dendograma  
    # Datos
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
    # Creamos uno para cada nivel del factor
    # Niveles a tomar
    Fac = glue("Data${factor}")
    Levs = eval(parse(text=Fac))
    # División
    T1 = "Data[,which(Levs == x)]" 
    
    # Martiz
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
            # Division por tejidos
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
  # Función: RX_DendoRelations
  #
  #   Crea una representación visual de las relaciones entre dos dendrogramas.
  #
  # Input:
  # -----
  #   - DendoData1 (objeto de clase "dendro"):
  #         Los datos del primer dendrograma.
  #   - DendoData2 (objeto de clase "dendro"):
  #         Los datos del segundo dendrograma.
  #   - size (vector numérico, por defecto c(0.65, 10, 0)):
  #         El tamaño relativo de las diferentes partes del gráfico. El vector
  #         debe tener tres elementos que representen el tamaño del panel 
  #         superior, el tamaño del panel de las líneas y el tamaño del panel 
  #         inferior, respectivamente.
  #
  # Output:
  # -----
  #   - (objeto de la clase "ggplot")
  #
  # Dependencias:
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
  # Función: RX_contrast
  #
  #   Esta función realiza un análisis de expresión diferencial utilizando 
  #   el paquete "limma". 
  #   Calcula los resultados específicos para un conjunto de datos de expresión
  #   "expmatrix", la matriz de diseño "design" y el contraste específico "C".
  #   Finalmente, ajusta los p-valores utilizando el método "BH" y
  #   devuelve los resultados.
  #
  # Input:
  # -----
  #   - expmatrix (data.frame o matriz):
  #         El conjunto de datos de expresión utilizado en
  #         el análisis de contraste.
  #   - design (matriz):
  #         Matriz de diseño experimental utilizado en el análisis de contraste.
  #   - C (character):
  #         Contraste específico a computar.
  #   - seed (valor numérico, por defecto 1808):
  #         La semilla utilizada para generar números aleatorios en el análisis.
  #
  # Output:
  # -----
  #   - (objeto de la clase "MArrayLM")
  #       Resultados del análisis de contraste.
  #
  # Dependencias:
  # -----
  #   - library(limma)
  #
  # ****************************************************** 
  
  set.seed(seed)
  # Matriz de contrastes
  contMatrix <- makeContrasts(contrasts = C, levels = design)
  # Ajusta un modelo lineal
  fit <- lmFit(expmatrix, design)
  fit2 <- contrasts.fit(fit, contMatrix)
  # Realiza el análisis de expresión diferencial
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
  # Función: RX_DiffExp
  #
  #   Esta función realiza un análisis de expresión diferencial utilizando
  #   el paquete "limma". 
  #   Esta función toma como entrada un objeto phenoData que contiene 
  #   información fanotipica de las muestras y un objeto expressData que 
  #   contiene los datos de expresión.
  #   Calcula los contrastes especificados en "C" para los datos de 
  #   expresión "expressData" y los datos fenotípicos "phenoData".
  #   El análisis puede realizarse
  #   para estudios de tipo "Array" o "RNA-seq", y se pueden incluir covariables para el diseño.
  #   Devuelve una lista de objetos que contienen los resultados de cada contraste.
  #
  # Input:
  # -----
  #   - phenoData (data.frame):
  #         Los datos fenotípicos.
  #   - expressData (data.frame o matriz):
  #         Los datos de expresión.
  #   - C (character vector):
  #         Vector de los contrastes a computar.
  #   - var (carácter o vector de caracteres, por defecto "Group"):
  #         El nombre o los nombres de las variables en "phenoData" que se
  #         utilizarán para construir los contrastes.
  #         Si se proporcionan dos nombres, se creará una interacción entre
  #         las variables.
  #   - studyType (carácter, por defecto "Array"):
  #         El tipo de estudio, puede ser "Array" o "RNA-seq". si se trata de
  #         "RNA-seq", se aplicará el método "voom" a los datos de expresión.
  #   - covar (carácter, opcional):
  #         El nombre de la variable en "phenoData" que se utilizará como
  #         covariable en el diseño.
  #   - seed (valor numérico, por defecto 1808):
  #         La semilla utilizada para generar números aleatorios en el análisis.
  # Output:
  # -----
  #   - (list)
  #     Lista de resultados por contraste. Cada objeto es de la clase "MArrayLM"
  #     y contiene los resultados del análisis de expresión diferencial para un
  #     contraste específico.
  #
  # Dependencias:
  # -----
  #   - library(limma)
  #   - library(Biobase)
  #   - library(stringr)
  #
  # ******************************************************  
  set.seed(seed)
  
  if(length(var) == 2){
    # Crear un contraste de interacción entre las dos variables de phenoData
    contrast = interaction(phenoData[,var[1]], phenoData[,var[2]])
  }else{
    contrast = phenoData[,var]
  }
  # Evalur la presencia de covariables
  if(is.null(covar) | ! is.factor(phenoData[, covar])){
    # Crear una matriz de diseño sin covariables
    design = model.matrix(~0 + contrast)
    colnames(design) <- str_replace(levels(contrast) , "-", "_")
  }else{
    # Crear una matriz de diseño con el contraste y la covariable
    cov = factor(make.names(phenoData[, covar]))
    design = model.matrix(~ 0 + contrast + cov)
    colnames(design) <- c(str_replace(levels(contrast) , "-", "_"),
                          levels(cov)[-2])
  }
  # Si studyType es "RNA-seq", aplicar la transformación voom
  if (studyType == "RNA-seq"){
    expressData = voom(expressData, design = design, plot = F)
  }
  # Realizar el análisis de expresión diferencial para cada contraste
  Res <- sapply(C, function(x) RX_contrast(expressData, design, x),
                simplify = FALSE, USE.NAMES = TRUE)
  return(Res)
}


#------------- RX_DiffExpFinal
RX_DiffExpFinal <- function(Data,
                            var1="Group",
                            var2="Gender",
                            seed = 1808,
                            ...){
  # ******************************************************** 
  # Función: RX_DiffExpFinal
  #
  #   Esta función realiza un conjunto de análisis de expresión diferencial 
  #   para lso escenarios de un estudio utilizando el paquete "limma".
  #   Realiza análisis de expresión diferencial para diferentes tejidos en un
  #   conjunto de datos.Para cada uno de estos calcula los contrastes  
  #   especificados con la función "RX_DiffExp"  y combina los resultados en 
  #   una lista. 
  #
  # Input:
  # -----
  #   - Data (objeto ExpressionSet o SummarizedExperiment):
  #         El conjunto de datos que contiene los datos de expresión y
  #         los datos fenotípicos.
  #   - var1 (carácter, por defecto "Group"):
  #         El nombre de la variable en los datos fenotípicos que se
  #         utilizará para construir los
  #         contrastes.
  #   - var2 (carácter, por defecto "Gender"):
  #         El nombre de la variable en los datos fenotípicos que se utilizará
  #         para construir los contrastes.
  #   - seed (entero, por defecto 1808):
  #         La semilla utilizada para la generación de números aleatorios.
  #   - ... :
  #         Otros argumentos que se pasarán a la función "RX_DiffExp".
  #
  # Output:
  # -----
  #   - Lista de los resultados de expresión diferencial.
  #     Para cada tejido, estudio y contraste obtenemos un objeto de la
  #     clase MArrayLM
  #
  # Dependencias:
  # -----
  #   - library(SummarizedExperiment)
  #   - library(Biobase)
  #   - library(glue)
  #
  # ******************************************************** 
  set.seed(seed)
  # Determinamos la proveniencia de los datos de expresión
  if(class(Data) == "ExpressionSet"){
    phenoData = pData(Data)
    expressData = exprs(Data)
    studyType = "Array"
  }else{
    phenoData = colData(Data)
    studyType = "RNA-seq"
    expressData = assay(Data)      
  }
  # Tomamos las variables
  Group = str_replace(levels(phenoData[,var1]) , "-", "_")
  Gender = levels(phenoData[,var2])
  Tissue = levels(phenoData$Tissue)
  # Separamos por tejidos
  phenoDatas = sapply(Tissue, function(x)
    phenoData[which(phenoData$Tissue == x),],
                      simplify = FALSE, USE.NAMES = TRUE)
  expressDatas = sapply(Tissue, function(x)
    expressData[,which(phenoData$Tissue == x)],
                        simplify = FALSE, USE.NAMES = TRUE)
  
  V1 = str_replace(levels(phenoData[,var1]) , "-", "_")
  V2 = levels(phenoData[,var2]) 
  # Contrastes
  C1 = combn(rev(Group), m=2, FUN = paste,collapse=" - ", simplify = TRUE)
  
  C2 = sapply(Gender, 
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
                           colors = c("Significativo logFC > 0" = "#CB326D", 
                                      "Significativo logFC < 0" = "#00a0ab",
                                      "Significativo" = "gray24",
                                      "No significativo" = "grey"),
                           pval_col = "adj.P.Val",
                           logFC_col = "logFC",
                           point_size = 1,
                           point_shape = 8,
                           label = "ENTREZID",
                           legendTitle = "Significancia",
                           ylab = "- log10 (FDR)",
                           xlab = "logFC"){
  # ****************************************************** 
  # Función: RX_VolcanoPlot
  #
  #   Esta función crea un VolcanoPlot que muestra la relación entre el 
  #   logaritmo del fold change (logFC) y el logaritmo negativo del valor  
  #   p ajustado (FDR) de un conjunto de datos. Los puntos se colorearán según  
  #   su significancia y se pueden etiquetar con identificadores específicos.
  #
  # Input:
  # -----
  #   - DF (data.frame):
  #         El conjunto de datos que contiene las columnas necesarias para el 
  #         gráfico, incluyendo el logFC, el valor p ajustado, y posiblemente 
  #         identificadores para etiquetar los puntos.
  #   - logFC_lim (valor numérico, por defecto 1):
  #         El límite utilizado para definir los puntos significativos basados
  #         en el logFC.
  #   - padj_lim (valor numérico, por defecto 0.05):
  #         El límite utilizado para definir los puntos significativos basados 
  #         en el valor p ajustado.
  #   - FDRtransform (función, por defecto log10):
  #         La función utilizada para transformar el valor p ajustado en el eje 
  #         y del gráfico.
  #   - colors (vector de colores, por defecto 
  #             c("Significativo logFC > 0" = "#CB326D", 
  #               "Significativo logFC < 0" = "#00a0ab",
  #               "Significativo" = "gray24",
  #               "No significativo" = "grey")):
  #         Los colores utilizados para representar los diferentes grupos de
  #         puntos en el gráfico.
  #   - pval_col (nombre de columna, por defecto "adj.P.Val"):
  #         El nombre de la columna que contiene los valores p ajustados en
  #         el conjunto de datos.
  #   - logFC_col (nombre de columna, por defecto "logFC"):
  #         El nombre de la columna que contiene los valores de logFC en
  #         el conjunto de datos.
  #   - point_size (valor numérico, por defecto 1):
  #         El tamaño de los puntos en el gráfico.
  #   - point_shape (valor numérico, por defecto 8):
  #         La forma de los puntos en el gráfico.
  #   - label (nombre de columna, por defecto "ENTREZID"):
  #         El nombre de la columna que contiene los identificadores 
  #         utilizados para etiquetar los puntos.
  #   - legendTitle (texto, por defecto "Significancia"):
  #         El título de la leyenda del gráfico.
  #   - ylab (texto, por defecto "- log10 (FDR)"):
  #         Etiqueta del eje y del gráfico.
  #   - xlab (texto, por defecto "logFC"):
  #         Etiqueta del eje x del gráfico.
  #
  # Output:
  # -----
  #   - (objeto de la clase "ggplot"):
  #       Volcano obtenido con ggplot.
  #
  # Dependencias:
  # -----
  #   - library(ggplot2)
  #   - library(ggrepel)
  #
  # ****************************************************** 
  
  # Evaluamos la significancia
  DF$Sig =  factor(sapply(seq(nrow(DF)), 
                          function(i)
                            ifelse(DF[i,pval_col]>padj_lim, "No significativo",
                                   ifelse(DF[i,logFC_col] < -logFC_lim, 
                                          "Significativo logFC < 0",
                                          ifelse(DF[i,logFC_col] > logFC_lim,
                                                 "Significativo logFC > 0",
                                                 "Significativo" )))
                          ))
  
  # Preparamos los datos
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
                      label = ifelse(DF$Sig %in% c("Significativo logFC < 0",
                                                   "Significativo logFC > 0"), 
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
  # Función: RX_getSE
  #
  #   Esta función calcula el error estándar (SE) de los coeficientes estimados
  #   en un objeto de ajuste limma mediante dos métodos.
  #
  # Input:
  # -----
  #   - Data (objeto de la clase "MArrayLM"):
  #         El objeto de ajuste limma que contiene los coeficientes estimados
  #         y la información necesaria para calcular el SE.
  #
  # Output:
  # -----
  #   - (data.frame):
  #       Un data frame que contiene los coeficientes estimados y sus
  #       respectivos SE.
  #
  # Dependencias:
  # -----
  #   - library(limma)
  #
  # ******************************************************  
  
  # Método 1:
  print(class(Data))
  TopTab <- topTable(Data, number = "all", confint=TRUE,
                     adjust.method = "fdr", sort.by = "none")
  # Calculamos el SE
  TopTab[, "SE"] <- (TopTab[, "CI.R"] - TopTab[, "CI.L"])/ 3.92
  
  # Método 2:
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
  # Función: RX_MetaAnalysis_prep
  #
  #   Esta función prepara los datos para realizar un metaanálisis de expresión 
  #   génica utilizando los resultados de varios estudios. Se extraen los 
  #   datos de los tejidos y contrastes indicados,se calcula el error estándar 
  #   (SE) correspondiente a los coeficientes estimados, y se crea una matriz 
  #   de logFC y una matriz de SE para su uso en el metaanálisis.
  #
  # Input:
  # -----
  #   - Studies (vector de caracteres):
  #         Los nombres de los estudios incluidos en el metaanálisis.
  #   - Tissue (vector de caracteres, uno o varios de c("SAT", "VAT", "LSAT")):
  #         Los tejidos a incluir en el metaanálisis. Puede ser uno o varios tejidos.
  #   - Contrast (vector de caracteres):
  #         Los contrastes a utilizar en el metaanálisis.
  #   - Datas (lista):
  #         Los datos de expresión génica de cada estudio y tejido,
  #         organizados en una lista.
  #   - SE_coef (vector de caracteres, uno de c("SE", "coefSE")):
  #         Método de obtención del error estándar (SE) en los datos de cada 
  #         estudio.
  # Output:
  # -----
  #   - (lista):
  #       Una lista que contiene la matriz de logFC, la matriz de SE y la
  #       información de anotación para el metaanálisis.
  #
  # Dependencias:
  # -----
  #   - library(dplyr)
  #   - library(limma)
  #
  # ****************************************************** 
  
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


