
setMethod("show", "blum", function(object){
  cat("An object of class blum with",nrow(fData(object)),"analytes:","\n")
  cat("\t", as.character(head(fData(object)$analyte, 3)),"...", as.character(tail(fData(object)$analyte, 3)),"\n")
  cat(nrow(exprs(object)), "measures of expression in", nrow(pData(object)),"wells, on", length(unique(pData(object)$plate)), "plates.","\n")
  cat("And slots:", names(getSlots("blum")),"\n")
})


setMethod("head", signature=c("blum"), function(x){
  cat("@featureData:\n")
  print(head(fData(x)))
  cat("@phenoData:\n")
  print(head(pData(x)))
  cat("@exprs:\n")
  print(exprs(x))
})


setMethod("head", signature=c("slum"), function(x){
  cat("@formula:\n")
  print(formula(x))
  cat("@featureData:\n")
  print(head(fData(x)))
  cat("@fit:\n")
  print(head(fit(x)))
  cat("@phenoData:\n")
  print(head(pData(x)))
  cat("@exprs:\n")
  print(exprs(x)[1:5, 1:4])
})


setMethod("pData", "blum", function(object){
  return(pData(object@phenoData))
})
setReplaceMethod("pData", signature("blum", "data.frame"), function(object, value){
  object@phenoData@data<-value
  return(object)})

setMethod("fData", "blum", function(object){
  return(pData(object@featureData))
})
setReplaceMethod("fData", signature("blum", "data.frame"), function(object, value){
  object@featureData@data<-value
  return(object)})


setMethod("exprs", "blum", function(object){
  return(object@exprs)
})
setReplaceMethod("exprs", signature("blum", "data.table"), function(object, value){
  object@exprs<-value
  return(object)})


setGeneric("fit", function(object, ...) standardGeneric("fit"))
setMethod("fit", "slum",function(object){
  return(object@fit)
})

#' Concentration accessor
#' 
#' Concentration matrix accessor for slum objects.
#' 
#' @param object A \code{slum} object.
#' 
#' @details
#' The function access the concentration matrix in the \code{assayData} slot of
#' a \code{slum} object
#' 
#' @seealso \code{\link{slum}}
#' @author Renan Sauteraud
#' 
#' @aliases
#' concentration
#' concentration,slum-method
#' @exportMethod concentration
setGeneric("concentration" ,function(object) standardGeneric("concentration"))
setMethod("concentration", "slum", function(object){
  #return(assayData(object)$concentration)
  return(assayDataElement(object, "concentration"))
})


setGeneric("subset") #create new S4 generic
setMethod(f="subset", signature="blum", definition=function(x, subset, select, ...){
  if(missing(subset)){
    fdata_call <- substitute(TRUE)
  } else{
    fdata_call <- substitute(subset)
  }
  if(missing(select)){
    pdata_call <- substitute(TRUE)
  } else{
    pdata_call <- substitute(select)
  }
  fdata_rows <- eval(fdata_call, fData(x), parent.frame())
  pdata_rows <- eval(pdata_call, pData(x), parent.frame())
  fData(x) <- fData(x)[fdata_rows,]
  pData(x) <- pData(x)[pdata_rows,]
  #exprs should be subset only based on sample_id and analyte
  #exprs(x) <- exprs(x)[eval(fdata_call, exprs(x), parent.frame()) & eval(pdata_call, exprs(x), parent.frame()),]
  exprs(x) <- exprs(x)[analyte %in% fData(x)[, "analyte"] & sample_id %in% pData(x)[, "sample_id"], ]
  x
})

#' @export melt
#' @export
melt.slum <- function(data, ..., na.rm, value.name) {
  mslum <- data.table(melt(exprs(data)))
  setnames(mslum, colnames(mslum), c("analyte","sample_id",tolower(data@unit)))
  mslum <- merge(mslum, pData(data), by="sample_id")
  mslum <- merge(mslum, fData(data), by="analyte")
  return(mslum)
}

#' @export
melt.blum <- function(data, ..., na.rm, value.name) {
  mblum <- exprs(data)
  mblum <- merge(mblum, pData(data), by="sample_id")
  mblum <- merge(mblum, fData(data), by="analyte")
  return(mblum)
}




#--------
# formula<-
# setter
setGeneric("formula<-", function(object, value) standardGeneric("formula<-"))
setReplaceMethod("formula", signature("slum", "formula"), function(object, value){
  object@formula<-value
  return(object)})


#' Plot the standard curves
#' 
#' Create a geom_line object to add the standard curves to a plot.
#' 
#' @param object A \code{slum}. The object that contains the experiment to plot.
#' @param n A \code{numeric}. The number of points used to draw the curves.
#' @param na.rm A \code{logical}. Set to TRUE if some of the calculated 
#'  concentrations are NA.
#' @param mapping,data,stat,position,... Additional parameters to be passed 
#'  to geom_line
#' 
#' @details
#' See ggplot2's documentation for more information regarding the additional 
#' arguments.
#'   
#' @seealso \code{\link{geom_line}}
#' @author Renan Sauteraud
#'  
#' @aliases
#' geom_sc,slum-method
#' @importFrom ggplot2 geom_line
#' @exportMethod geom_sc
setGeneric("geom_sc", function(object, n=100, data = NULL, mapping = aes(x=concentration, y=mfi),
				stat = "identity", position = "identity", 
				na.rm = FALSE, ...) standardGeneric("geom_sc"))

# setMethod("geom_sc", "slum",
#           function(object, n=100, mapping = aes(x=concentration, y=mfi),
#                    stat = "identity", position = "identity", 
#                    na.rm = FALSE, ...){            
#             fit<-object@fit
#             fit$plate <- as.character(fit$plate)
#             # Extract the coefficients per plate/analyte
#             # Split by plate then analyte                                
#             sfit<-lapply(split(fit,fit$plate),function(x){split(x,x$analyte)})
#             df.sc<-lapply(sfit,lapply,.df_sc_fit,object@formula, n, na.rm)
#             df.sc<-as.data.frame(do.call("rbind",lapply(df.sc,function(x)do.call("rbind",x))))
#             
#             ret<-geom_line(data=df.sc, mapping=mapping,
#                            na.rm=na.rm, ...)
#             return(ret)
#           })

setMethod("geom_sc", "slum",
  function(object, n = 100, mapping = aes(x = concentration, y = mfi), 
           stat = "identity", position = "identity", na.rm = TRUE, ...){
    dt_fit <- data.table(object@fit)
    dt_conc <- dt_fit[, exp(seq(0, log(max(concentration, na.rm = TRUE)), length.out =100)), by = c("plate", "analyte")]
    dt_conc <- merge(dt_conc, unique(dt_fit[, list(plate, analyte, b,c,d,e,f)]), by = c("plate", "analyte"))
    setnames(dt_conc, "V1", "x")
    dt_conc[, mfi := eval(parse(text=object@formula[3]))]
    dt_conc[, mfi := exp(mfi)]
    dt_conc[, c(letters[2:6]) := NULL]
    setnames(dt_conc, "x", "concentration")
    ret<-geom_line(mapping=mapping, data=dt_conc, stat = stat, 
                   position = position, na.rm=na.rm, ...)
  })

# .df_sc_fit<-function(df, formula, n=100, na.rm=TRUE){
#   # Concentration from min to max
#   x<-exp(seq(0,log(max(df$concentration, na.rm = na.rm)),length.out=n))
#   coef<-as.numeric(df[1,c('b','c','d','e','f')])
#   b<-coef[1];c<-coef[2];d<-coef[3];e<-coef[4];f<-coef[5];
#   # Evalulate formula
#   mfi<-eval(parse(text=formula[3]))
#   # Check whether the formula is fit on the log scale
#   if(substr(as.character(formula[2]),1,3)=="log")
#     mfi <- exp(mfi)  
#   # basic dataframes with plate, filename, well, concentration, mfi
#   newdf<-data.frame(plate=rep(df$plate[1],n),analyte=rep(df$analyte[1],n),concentration=x,mfi=mfi)
#  return(newdf)
# }

#' Plot a plate layout
#' 
#' Visualisation of a plate's sample. This draws a plot of a selected plate 
#' where each graphical element can be changed to display sample information.
#' 
#' @param object A \code{blum} or \code{slum}. The object containing the 
#'  experiment to plot.
#' @param plate_name A \code{character}. The name of the plate to plot.
#' @param ... Additional arguments to be passed to ggplot2's \code{geom_polygon}
#'  function.
#'  
#' @details
#' The information passed in the additional aguments must be available in the
#' \code{phenoData} slot of the object.
#' 
#' @examples
#' #plot_layout(blum, "plate1", fill="sample_type", col = "concentration")
#'  
#' @seealso \code{\link{geom_polygon}}
#' @author Renan Sauteraud
#' 
#' @aliases
#' plot_layout,blumORslum-method
#' @importFrom ggplot2 ggplot geom_polygon aes aes_string labs 
#' @importFrom ggplot2 theme element_blank facet_grid
#' @importFrom grid unit
#' @exportMethod plot_layout
setGeneric("plot_layout", function(object, plate_name=NULL, ... ) standardGeneric("plot_layout"))
setMethod("plot_layout", "blumORslum", function(object, plate_name=NULL, ...){
  pd<-pData(object)
  plateNames<-levels(pd$plate)
  if(is.null(plate_name)){
    plate_name<-levels(pd$plate)[1]
    warning("Plate name not specified '",plate_name,"' will be displayed")
  } else if(!plate_name%in%plateNames) {
    stop("'",plate_name,"' is not a valid plate name. Available plate names for this object are: ",paste(plateNames, collapse=", "))
  }
  pd<-subset(pd,plate==plate_name)
  df<-data.frame(well2coords(pd$well), pd)
  df<-cbind(df, x=rep(1, nrow(df)), y=rep(1, nrow(df)))
  
  df2<-lapply(df$well,function(well){
    angle <- seq(-pi, pi, length = 50);
    xx = sin(angle); yy = cos(angle);
    data.frame(well,xx,yy)
    })
  df2<-do.call("rbind",df2)
  
  df.combine<-merge(df,df2,by="well")
  p<-ggplot(df.combine, aes_string("x", "y")) + 
    geom_polygon(mapping=aes_string(x="xx", y="yy", ...),data=df.combine)
  # geom_point(aes(x=x,y=y, color=sample_type))+theme_bw()
  p<-p+labs(title=plate_name)+theme(line=element_blank(), axis.text=element_blank(), axis.title=element_blank(), panel.margin = unit(0, "lines"))+facet_grid(row~col)
  return(p)
})


well2coords<-function(well_id){
  row<-substr(well_id, 1,1)
  col<-as.numeric(substr(well_id, 2,3))
  return(list(row=row, col=col))
}

setGeneric("getCoeffs", function(object, plate=NULL, analyte=NULL) standardGeneric("getCoeffs"))
setMethod("getCoeffs", "slum", function(object, plate=NULL, analyte=NULL){
  if(is.null(plate) | is.null(analyte)){
    stop("Missing argument  'plate' or 'analyte'")
  }
  df<-unique(object@fit[,c("plate","analyte","b","c","d","e","f")])
  df<-df[df$plate==plate & df$analyte==analyte, c("b","c","d","e","f")]
  if(nrow(df)==0){
    stop("No match found for the given 'plate' or 'analyte'.")
  }
  return(df)
  #return(as.numeric(df))
})

setMethod("merge", signature=c("blum", "blum"), function(x, y, ...){#obj1, obj2, ...){
  if(missing(y)){
    return(x)
  }
  list_obj<-c(list(x, y), list(...))
  if(.is.mergeable(list_obj)){
    phenoData<-as(do.call("rbind", lapply(list_obj, pData)), "AnnotatedDataFrame")
    featureData<-x@featureData
    exprs<-do.call("rbind", lapply(list_obj, exprs)) # Mighy be tricky if colnames are =/=
    ret<-new("blum", phenoData=phenoData, featureData=featureData, exprs=exprs)
    return(ret)
  }
})

.is.mergeable<-function(list_obj){
  # Plate names must be unique in the merged object
  allPlates<-unlist(lapply(list_obj, function(x){ unique(pData(x)$plate) }))
  if(length(unique(allPlates))!=length(allPlates)){
    duplates<-allPlates[duplicated(allPlates)]
    stop("Some plates are not unique. Merging an object with itself?")
  }
  # fData are the same 
  fds<-lapply(list_obj, function(x){ fData(x) })
  fds<-lapply(fds, function(x){ x[with(x, order(bid)),] }) # sort by bid
  fdMatch<-all(unlist(lapply(fds, identical, fds[[1]]))) # compare all df to the first one
  if(fdMatch==FALSE){
    stop("The analytes/bid of the objects to merge are different")
  }
  # layout is identical
  pds<-lapply(list_obj, function(x){ pData(x) })
  pds<-lapply(pds, function(x){ x[,!colnames(x)%in%c("plate", "filename", "sample_id", "center")] })
  pds<-lapply(pds, function(x){ x[with(x, order(well)),] })
  pdMatch<-all(unlist(lapply(pds, identical, pds[[1]])))
  if(pdMatch==FALSE){
    stop("The phenotype data of the objects to merge are different. Different layout?")
  }
  
  return(TRUE)
}

#' set_center
#' 
#' Add information regarding which center or lab produced the data. This is 
#' useful when multiple source ran the same experiment.
#' 
#' @param object A \code{blum} or \code{slum} object.
#' @param center_name A \code{character}. The name of the center or lab that 
#'  generated the data in the object.
#'  
#' @seealso \code{blum}, \code{slum}
#' @author Renan Sauteraud
#' 
#' @aliases
#' set_center,blum,character-method
#' set_center,slum,character-method
#' @import data.table
#' @exportMethod set_center
setGeneric("set_center", function(object, center_name) standardGeneric("set_center"))
setMethod("set_center", signature=c("blum", "character"), function(object, center_name){
  pd<-pData(object)
  dt<-data.table(exprs(object))
  if("center"%in%colnames(pd)){
    pd$plate<-gsub(paste0(pd$center[1],"_"), "",  pd$plate)
    dt[,plate:=gsub(paste0(center[1],"_"), "", plate)]
  }
  pd$center<-rep(center_name, nrow(pd))
  pd$plate<-paste(center_name, pd$plate, sep="_")
  pd$sample_id<-paste(sapply(strsplit(as.character(pd$filename), split="\\."), function(x){x[-length(x)]}),
        pd$plate, pd$well, sep="_")

  dt[,center:=rep(center_name, nrow(dt))]
  dt[,plate:=paste(center_name, plate, sep="_")]
  dt[,sample_id:=paste(sapply(strsplit(filename, split="\\."), function(x){x[-length(x)]}), plate, well, sep="_")]
  pData(object)<-pd
  exprs(object)<-dt
  return(object)
})

setMethod("set_center", signature=c("slum", "character"), function(object, center_name){
  pd<-pData(object)
  cNames<-colnames(exprs(object))
  if("center"%in%colnames(pd)){
    cNames<-gsub(pd$center[1], center_name, cNames)
    pd$plate<-gsub(pd$center[1], center_name,  pd$plate)
    pd$sample_id<-paste(sapply(strsplit(as.character(pd$filename), split="\\."), function(x){x[-length(x)]}),
        pd$plate, pd$well, sep="_")
  } else{
    for(plate in unique(pd$plate)){
      cNames<-gsub(paste0(plate,"_"), paste(center_name, "_", plate, "_",  sep=""), cNames)
    }
    pd$plate<-paste(center_name, pd$plate, sep="_")
    pd$sample_id<-paste(sapply(strsplit(as.character(pd$filename), split="\\."), function(x){x[-length(x)]}),
        pd$plate, pd$well, sep="_")
  }
  colnames(exprs(object))<-cNames
  colnames(assayData(object)$concentration)<-cNames
  pd$center<-rep(center_name, nrow(pd))
  pData(object)<-pd
  return(object)
})


