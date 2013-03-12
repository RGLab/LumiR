# show method:
#   Shows class, slots, number of analytes, total number of measures
#
setMethod("show", "blum", function(object){
  cat("An object of class blum with",nrow(fData(object)),"analytes:","\n")
  cat("\t", as.character(head(fData(object)$analyte, 3)),"...", as.character(tail(fData(object)$analyte, 3)),"\n")
  cat(nrow(exprs(object)), "measures of expression in", nrow(pData(object)),"wells, on", length(unique(pData(object)$plate)), "plates.","\n")
  cat("And slots:", names(getSlots("blum")),"\n")
})

# pData method:
#   phenoData accessor
setMethod("pData", "blum", function(object){
  return(pData(object@phenoData))
})
setReplaceMethod("pData", signature("blum", "data.frame"), function(object, value){
  object@phenoData@data<-value
  return(object)})

## contains information about analytes
setMethod("fData", "blum", function(object){
  return(pData(object@featureData))
})
setReplaceMethod("fData", signature("blum", "data.frame"), function(object, value){
  object@featureData@data<-value
  return(object)})

# exprs accessor for bead level data a la eSet
setMethod("exprs", "blum", function(object){
  return(object@exprs)
})


# fit accessor for standard curve fitting information
setGeneric("fit", function(object, ...) standardGeneric("fit"))
setMethod("fit", "slum",function(object){
  return(object@fit)
})
setGeneric("concentration" ,function(object,...) standardGeneric("concentration"))
setMethod("concentration", "slum", function(object){
  return(assayData(object)$concentration)
})

# Subset method to subset a la eSet
setMethod("[","blum",
          function(x,i,j,..., drop=FALSE)
          {
            if(!missing(i))
            {
              #Subset the samples
              bdata<-exprs(x)[j]
              #Subset the analytes
              bdata<-lapply(exprs(x),"[",i)              
            }
            else
            {
              #Subset the analytes         
              bdata<-lapply(exprs(x),"[",i)              
            }            
            newSet<-new('blum'
                        ,exprs=bdata
                        ,phenoData=x@phenoData[j,]
                        ,featureData=x@featureData[i,])
            newSet            
          })

# I add the pheno and feature information by default we could add an option for this
# Need to sort this out with reshape2
# Also needs some work
setGeneric("melt",function(x,...){
  standardGeneric("melt")
})


setMethod("melt","blum",
function(x)
  {
  # Use the melt function in reshape2
  df<-reshape2::melt(exprs(x))
  names(df)<-c("FL","analyte","sample_id")
  df<-merge(df,pData(x),by=c("sample_id"))
  df<-merge(df,fData(x),by="analyte")
  return(df)
  
})

setMethod("melt","slum",
          function(x)
          {
            # Use the melt function in reshape2
            df<-reshape2::melt(exprs(x))
            names(df)<-c("analyte","sample_id",tolower(x@unit))            
            ## merge all information
            df<-merge(df,pData(x),by=c("sample_id"))
            df<-merge(df,fData(x),by="analyte")
            return(df)
          })

#--------
# formula
# getter
#setGeneric("formula", function(object, ...) standardGeneric("formula"))
#setMethod("formula", "slum", function(object){return(object@formula)})
#--------
# formula<-
# setter
setGeneric("formula<-", function(object, value, ...) standardGeneric("formula<-"))
setReplaceMethod("formula", "slum", function(object, value){object@formula<-value; object})


setGeneric("geom_sc", function(object, n=100, data = NULL, mapping = aes(x=concentration, y=mfi),
				stat = "identity", position = "identity", 
				na.rm = FALSE, ...) standardGeneric("geom_sc"))

setMethod("geom_sc", "slum",
          function(object, n=100, mapping = aes(x=concentration, y=mfi),
                   stat = "identity", position = "identity", 
                   na.rm = FALSE, ...)
          {            
            
            fit<-object@fit
            # Extract the coefficients per plate/analyte
            # Split by plate then analyte                                
            sfit<-lapply(split(fit,fit$plate),function(x){split(x,x$analyte)})
            df.sc<-lapply(sfit,lapply,.df_sc_fit,object@formula)
            df.sc<-as.data.frame(do.call("rbind",lapply(df.sc,function(x)do.call("rbind",x))))
            
            ret<-geom_line(data=df.sc, mapping=mapping,
                           na.rm=na.rm, ...)
            return(ret)
          })

.df_sc_fit<-function(df,formula,n=100)
{
  # Concentration from min to max
  x<-exp(seq(0,log(max(df$concentration)),length.out=n))
  coef<-as.numeric(df[1,c('b','c','d','e','f')])
  # Should be able to parse this automatically
  b<-coef[1];c<-coef[2];d<-coef[3];e<-coef[4];f<-coef[5];
  # Evalulate formula
  mfi<-eval(parse(text=formula[3]))
  # Check whether the formula is fit on the log scale
  if(substr(as.character(formula[2]),1,3)=="log")
    mfi <- exp(mfi)  
  # basic dataframes with plate, filename, well, concentration, mfi
  newdf<-data.frame(plate=rep(df$plate[1],n),analyte=rep(df$analyte[1],n),concentration=x,mfi=mfi)
 return(newdf)
}

setGeneric("plot_layout", function(object, plate_name=NULL,... ) standardGeneric("plot_layout"))

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
  p<-ggplot(df.combine, aes(x, y)) + 
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

