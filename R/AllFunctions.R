#' Create templates
#' 
#' Create files to be used as mapping files. The function will guess as much 
#' information as possible based on the data files. The created files may be
#' incomplete and require user modifications but they ensure that the required
#' columns are present.
#' 
#' @param path A \code{character}. The pathname of an experiment folder. This is
#' the same path as the one that will be used in \code{read.experiment}.
#' @param templates A \code{character} vector. Can be any combination of 
#' "layout", "analyte" and "phenotype". A mapping file will be created for each
#' selected template.
#' @param write A \code{logical}. If TRUE, write the template files in the 
#'  given \code{path}. When \code{setup_templates} is called through 
#'  \code{read.experiment}, write is set to FALSE and the template files are not
#'  saved.
#' 
#' @details
#' If a mapping file specified by the templates list already exists, the function
#' will not overwrite the existing file. The files should be removed before 
#' running the function in order to create new mapping files.
#' 
#' @return An invisible \code{list} of \code{data.frame} to be used internally by
#' \code{read.experiment} if one or more mapping files are missing from the given
#' path.
#' 
#' @author Renan Sauteraud
#' 
#' @seealso \code{\link{read.experiment}}
#' 
#' @export
#' @importFrom flowCore read.FCS
#' @importMethodsFrom flowCore exprs
#' @importFrom XML xmlRoot xmlTreeParse xmlAttrs xmlValue xmlSApply xmlSize xmlName
#' @importFrom tools list_files_with_exts
#' 
setup_templates<-function(path, templates=c("layout", "analyte", "phenotype"), write=TRUE){
  dfList <- list()
  analyte.file <- list.files(path, pattern = "analyte", full.names = TRUE)
  layout.file  <- list.files(path, pattern = "layout", full.names = TRUE)
  pheno.file   <- list.files(path, pattern = "phenotype", full.names = TRUE)
  plates <- list.dirs(path, recursive = FALSE)

  if(length(list_files_with_exts(plates[1], exts="lxb")>0)){
    type<-"LXB"; typeExt<-"lxb";
  } else if(length(list_files_with_exts(plates[1], exts="xml")>0)){
    type<-"BIOPLEX"; typeExt<-"xml";
  } else {
    type<-"XPONENT"; typeExt<-"csv"
  }
  
  all.files<-unlist(lapply(plates,list_files_with_exts,exts=typeExt))  

  #pData
  if("phenotype"%in%templates){
    if(length(pheno.file)>0){
      warning("The phenotype mapping file already exists, remove it to setup a template for it")
    } else {
      if(type=="BIOPLEX"){
        wells<-.getBioplexWellsID(all.files)
        wellsPerFile<-sapply(all.files, function(x){
          xmlSize(xmlRoot(xmlTreeParse(x))[["Wells"]])
        })
        #fNames<-rep(unlist(lapply(all.files, lapply, function(x)tail(strsplit(x,"/")[[1]],1))), wellsPerFile)
        fNames <- rep(basename(all.files), wellsPerFile)
        plate<-rep(plates, wellsPerFile)
        plate <- basename(plate)
      } else {
        len<-lapply(plates, function(x){length(list_files_with_exts(x,exts=typeExt))})
        fNames <- basename(all.files)
        plate<-rep(plates, len)
        plate <- basename(plate)
        if(length(plate) > length(fNames)){
          warning("Less filenames than expected. A plate folder may be missing a data file.")
        }
        if(type=="XPONENT"){
          wells<-.getXponentWellsID(fNames)
        } else if(type=="LXB"){
          wells<-.getLXBWellsID(fNames)
        } else {
          wells<-.getBioplexWellsID(fNames)
        }
      }
      noExt<-gsub(paste0(".",typeExt), "", fNames)
      sample_id<-paste(noExt, plate, wells, sep="_")
      phenotype.df<-data.frame(plate=plate,filename=fNames,well=wells, sample_id=sample_id, row.names=sample_id)
      dfList[["pheno"]]<-phenotype.df
      if(write){
        write.csv(phenotype.df, file=paste0(path,"phenotype.csv"), row.names=FALSE)
      }
    }
  }
    
  #layout
  if("layout"%in%templates){
    if(length(layout.file)>0){
      warning("The layout mapping file already exists, remove it to setup a template for it")
    } else {
      if(type == "BIOPLEX"){
        layout.df <- .getBioplexLayout(all.files[[1]]) #Because all layout have to be the same
        dfList[["layout"]]<-layout.df
      } else {
        wells<-paste0(LETTERS[1:8], rep(seq(1,12), each=8))
        sample_type<-rep("unknown", length(wells))
        concentration<-rep(NA, length(wells))
        layout.df<-data.frame(well=wells, sample_type=sample_type, concentration=concentration)
        dfList[["layout"]]<-layout.df
      }
      if(write){
        write.csv(layout.df, file=paste0(path,"layout.csv"), row.names=FALSE)
      }
    }
  }

  #analyte
  if("analyte"%in%templates){
    if(length(analyte.file)>0){
      warning("The analyte mapping file already exists, remove it to setup a template for it")
    } else {
      if(type=="LXB" & length(list_files_with_exts(path, exts="lxd")>0)){
          lxdFile<-list_files_with_exts(path, exts="lxd")
          analyte.df<-.read.lxd(lxdFile[1])
      } else if(type=="BIOPLEX"){
        analyte.df <- .getBioplexAnalytes(all.files[1])
      } else{
        if(type=="XPONENT"){
          BIDs<-.getXponentBID(all.files[1])
        } else if(type=="LXB"){ 
          BIDs<-.getLXBBID(all.files[1])
        }
        BIDs<-BIDs[BIDs!=0]
        analyte<-paste0(rep("unknown", length(BIDs)), BIDs)
        analyte.df<-data.frame(bid=BIDs, analyte=analyte)
      }
        dfList[["analyte"]]<-analyte.df
      if(write){
        write.csv(analyte.df, file=paste0(path,"analyte.csv"), row.names=FALSE)
      }
    
    }
  }
  return(invisible(dfList))
}

#.getLXDPheno <- function(lxdFiles){
#  #1lxd/plate or just 1lxd
#  for(lxdF in lxdFiles){
#    xml <- xmlParse(lxdF)
#    locName <- sapply(root[["Plate"]][which(names(root[["Plate"]]) == "Well")], function(x){ xmlValue(x[["LocName"]])})
#  }

.getXponentBID<-function(firstFile){
  #sLine<-grep("[Ee]vent[Nn]o", readLines(firstFile[1], n=5))-1
  #con<-read.csv(firstFile, skip=sLine, header=TRUE);
  #BIDs<-sort(as.numeric(unique(con[,2])))
  dt<-fread(firstFile)
  BIDs<-unique(sort(as.numeric(dt[[colnames(dt)[2]]])))
  return(BIDs)
}
.getLXBBID<-function(firstFile){
  suppressWarnings(lxb<-read.FCS(firstFile))
  BIDs<-sort(as.numeric(unique(exprs(lxb)[,1])))
  return(BIDs)
}
.getBioplexBID<-function(firstFile){
  xml<-xmlTreeParse(firstFile)
  root<-xmlRoot(xml)
  ROIsNode <- root[["Wells"]][[1]][["RunSettings"]][["RegionsOfInterest"]]
  BIDs<-as.numeric(xmlSApply(ROIsNode, xmlAttrs))
  #BIDs<-as.numeric(xmlSApply(root[["Wells"]][[1]][["RunSettings"]][["RegionsOfInterest"]], xmlAttrs))
  return(BIDs)
}

.getBioplexAnalytes<-function(firstFile){
  xml<-xmlTreeParse(firstFile)
  root<-xmlRoot(xml)
  AnalytesNode <- root[["Samples"]][[1]][[1]][["Analytes"]]
  BIDs<-as.numeric(xmlSApply(AnalytesNode, xmlAttrs))
  analytes<-xmlSApply(AnalytesNode, function(x){xmlValue(x[["AnalyteName"]])})
  #BIDs<-as.numeric(xmlSApply(root[["Samples"]][[1]][[1]][["Analytes"]], xmlAttrs))
  #analytes<-xmlSApply(root[["Samples"]][[1]][[1]][["Analytes"]], function(x){xmlValue(x[["AnalyteName"]])})
  analyte.df<-data.frame(analyte=analytes, bid=BIDs)
  return(analyte.df)
}

.read.lxd<-function(filename){#Workset>Setup>Region
  xml<-xmlTreeParse(filename)
  root<-xmlRoot(xml)
  setupNode<-root[["Setup"]]
  regionsIdx<-which(xmlSApply(setupNode, xmlName)=="Region")
  nRegion<-length(regionsIdx)
  ldf<-vector('list', nRegion)
  for(idx in 1:nRegion){
    ldf[[idx]]<-xmlAttrs(setupNode[[regionsIdx[[idx]]]])[c("name", "id")]
  }
  mat<-do.call(rbind, ldf)
  colnames(mat)<-c("analyte", "bid")
  df<-as.data.frame(mat)
  df[,"bid"]<-as.numeric(levels(df$bid))[df$bid]
  return(df)
  #featureData<-as(df, 'AnnotatedDataFrame')
  #return(featureData)
}

.getBioplexLayout <- function(fName){
  xml <- xmlTreeParse(fName)
  root <- xmlRoot(xml)
  node <- root[["Samples"]]
  #sample_type <- tolower(unlist(xmlSApply(node, names)))
  mw <- xmlSApply(node, function(x){ xmlSApply(x,
                        function(xxx){ length(xxx[["MemberWells"]])
                  })})
  rep <- unlist(mw, recursive = TRUE, use.names = FALSE)
  nm <- unlist(sapply(mw, names), use.names = FALSE)
  sample_type <- rep(nm, rep)
  exp_concs <- unlist(xmlSApply(node, function(x){ xmlSApply(x,
                                      function(xx){ as.numeric(xmlValue(xx[["Analytes"]][[1]][["ExpectedConc"]])) })}))
  rep_cnt <- unlist(xmlSApply(node, function(x){ xmlSApply(x, 
                                    function(xx){as.numeric(xmlValue(xx[["Analytes"]][[1]][["ReplicateCount"]])) }) }))
  concentration <- rep(exp_concs, rep_cnt)
  well <- unlist(xmlSApply(node, function(x){ xmlSApply(x,
                                 function(xx){ xmlSApply(xx[["MemberWells"]],
                                 function(xxx){ paste(LETTERS[as.numeric(xmlAttrs(xxx)["RowNo"])], xmlAttrs(xxx)["ColNo"], sep="")}
                )})}))
  layout.df <- data.frame(well = well, sample_type = sample_type, concentration = concentration)
  #layout.df <- layout.df[gtools::mixedorder(ldf$well),] #if I need to increase readability of layout.csv
  return(layout.df)
}


####################
##    RESULTS     ##
####################
results.conc.CSV<-function(object, file="./concentrations.csv"){
  mbs<-melt(object)
  concentration<-c()
  for(i in 1:nrow(mbs)){
    coefs<-getCoeffs(object, mbs[i,"plate"], mbs[i, "analyte"])
    concentration<-c(concentration,as.numeric(object@inv(mbs[i, "mfi"], coefs)))
  }
  toWrite<-cbind(mbs[,c("plate", "well", "analyte", "mfi")], concentration)
  write.csv(toWrite, file=file, row.names=FALSE)
  
  return(invisible(toWrite))
}

results.curves.CSV<-function(object, file="./curves.csv"){
  bsfo<-object@formula[3]
  fList<-c()
  bsfi<-unique(object@fit[,c("plate", "analyte", "b","c","d","e","f")])
  bsfi[,3:7]<-round(bsfi[,3:7], 4)
  for(i in 1:nrow(bsfi)){
    fList<-c(fList,gsub("c",bsfi[i,"c"],gsub("d",bsfi[i,"d"],gsub("e",bsfi[i,"e"],gsub("f",bsfi[i,"f"],bsfo)))))
  }
  toWrite<-cbind(bsfi[,c("plate", "analyte")], Formula=fList)
  maxConc<-max(pData(object)$concentration, na.rm=TRUE)
  minConc<-min(pData(object)[pData(object)$concentration>0,]$concentration, na.rm=TRUE)
  write.csv(toWrite, file=file, row.names=FALSE)
  return(invisible(toWrite))
}

##
#' Write a slum object as ImmPort MBAA results
#' 
#' Format a \code{slum} object into a table as specified by ImmPort template:
#' 'https://immport.niaid.nih.gov/example_submission_packages/MBAA_Results.xls'.
#' 
#' @param object A \code{slum} object. The summarized experiment to report.
#' @param outfile A \code{character}. The name of the output file without extension.
#' @param type A \code{character}. Can be "csv" or "xls". Determines the type of 
#' file to be written. This will also be appended to the filename given in outfile.
#' @param concentration_unit A \code{character}. The unit of the concentration
#' to be reported in the document.
#' 
#' @author Renan Sauteraud
#' 
#' @seealso \code{\link{slum}}
#' 
#' @importFrom tools file_ext
#' @importFrom xlsx write.xlsx
#' 
writeMBAA <- function(object, outfile="./MBAA_results", type="csv", concentration_unit="pg/mL"){
  outfile <- paste(outfile, type, sep=".")
  pd <- pData(object)
  SourceID <- SourceIDType <- AssayID <- AssayGroup <- AnalyteName <- c()
  MFI <- ConcentrationValue <- c()
  MFICoordinate <- c()
  for(i in 1:nrow(pd)){ #loop on samples
    sources <- rep(pd[i, "sample_id"], dim(object)[1])
    sources_type <- rep(as.character(pd[i, "sample_type"]), dim(object)[1])
    AID <- rep(pd[i, "plate"], dim(object)[1])
    AG <- rep(pd[i, "center"], dim(object)[1])
    analytes <- rownames(exprs(object))
    mfis <- exprs(object)[,pd[i, "sample_id"]]
    concs <- concentration(object)[,pd[i, "sample_id"]]
    mfics <- rep(as.character(pd[i, "well"]), dim(object)[1])
    #
    SourceID <- c(SourceID, sources)
    SourceIDType <- c(SourceIDType, sources_type)
    AssayID <- c(AssayID, AID)
    AssayGroup <- c(AssayGroup, AG)
    AnalyteName <- c(AnalyteName, analytes)
    MFI <- c(MFI, mfis)
    ConcentrationValue <- c(ConcentrationValue, concs)
    MFICoordinate <- c(MFICoordinate, mfics)
  }
  len <- length(AssayID)
  ConcentrationUnit <- rep(concentration_unit, len)
  OUT <- data.frame( `Source ID` = SourceID, `Source ID Type` = SourceIDType, 
                     `Assay ID` = AssayID, `Assay Group ID` = AssayGroup,
                     `Analyte Name` = AnalyteName, MFI = MFI,
                     `Concentration Value` = ConcentrationValue, `Concentration Unit` = ConcentrationUnit,
                     `MFI Coordinate` = MFICoordinate)
  if(type=="csv"){
    write.csv(OUT, file=outfile, row.names=FALSE)
  } else if(type=="xls"){
    write.xlsx(OUT, file=outfile, row.names=FALSE)
  } else {
    warning("writeMBAA only accepts 'csv' and 'xls' output. The file will be csv.")
    write.csv(OUT, file=outfile, row.names=FALSE)
  }
  #return(OUT)
}


####################
##  SUMMARY       ##
####################
.fivePL_formula <- as.formula("log(mfi) ~ c + (d - c)/(1 + exp(b * (log(x) - log(e))))^f")
.fivePL_inverse <- function(y, b,c,d,e,f){
  exp(log(((d - c)/(log(y) - c))^(1/f) - 1)/b + log(e))
}

##
#' Summarize a blum object
#' 
#' Summarize a blum object into a slum object. Fit the standard curves and 
#' calculate median fluorescence intensities (MFI) and the corresponding 
#' concentrations for each well.
#' 
#' @param from A \code{blum} object. The object to summarize. The object must 
#' have information about the standards location and concentration in its 
#' \code{phenodata} slot. See the user guide for more information on the 
#' requirements.
#' @param type A \code{character}. "MFI"
#' 
#' @return An object of class \code{slum}. With two expression matrices for MFI 
#' and concentrations.
#'  
#' @author Renan Sauteraud
#' 
#' @seealso \code{\link{blum}}, \code{\link{slum}}
#' @export
#' @importFrom drc drm LL.5
#' @importFrom reshape2 melt acast
#' 
slummarize<-function(from,type="MFI"){
  dt<-exprs(from)
  dt<-dt[,as.double(median(fl)), by="sample_id,analyte"]
  setnames(dt, c("sample_id", "analyte", type))
  setkey(dt, sample_id)
  mat<-matrix(dt[,MFI], ncol=length(unique(dt[,sample_id])), dimnames=list(unique(dt[,analyte]),unique(dt[,sample_id])))
  mfiSet <- new("slum", formula=.fivePL_formula, inv=.fivePL_inverse)
  exprs(mfiSet)<-mat
  pData(mfiSet)<-pData(from)
  fData(mfiSet)<-fData(from)
  mfiSet@unit="MFI"
  
  melt_slum <- melt(mfiSet)  
  melt_slum <- melt_slum[concentration!=0 & sample_type=="standard"]
  if(nrow(melt_slum)==0){
    stop("The object does not contain any standard. Edit layout.csv to include standards location and concentration.")
  }
  fit <- melt_slum
  setkeyv(fit, c("plate","analyte"))
  coefnames <- letters[2:6]
  #drm cares about the order of the data points
  #however, the fit should be of similar quality
  fit[,  eval(coefnames):= {res <- drm(log(mfi) ~ concentration, fct=LL.5())$coefficients; as.list(res)}, by="plate,analyte"]
  fit[, calc_conc:=.fivePL_inverse(mfi, b,c,d,e,f)]
  fit[, p100_recov:=calc_conc/concentration*100]
  mfiSet@fit <- data.frame(fit)
  conc_mat <- .get_concentration_matrix(mfiSet)
  assayData(mfiSet) <- list(exprs=mat, concentration=conc_mat)
  return(mfiSet)
}

.get_concentration_matrix <- function(object){
  fit_dt <- data.table(object@fit)
  exprs_dt <- data.table(melt(exprs(object)))
  pd_dt <- data.table(pData(object))
  setnames(exprs_dt, c("analyte", "sample_id", object@unit))
  setkeyv(fit_dt, c("analyte", "sample_id"))
  setkeyv(exprs_dt, c("analyte", "sample_id"))
  setkeyv(pd_dt, c("sample_id"))
  bdt <- merge(pd_dt, exprs_dt, by="sample_id")[,c("analyte", "sample_id", "plate", object@unit), with=FALSE]
  bdt <- merge(bdt, unique(fit_dt[, list(analyte, plate, b, c, d, e, f)]), by=c("analyte","plate"))
  bdt <- bdt[, concentration:=object@inv(MFI, b, c, d, e, f)]
  conc_mat <- acast(bdt, analyte ~ sample_id, value.var="concentration")
} 

#' Plot Standard curves
#' 
#' Wrapper around the geom_sc method for people who don't want to melt and 
#' ggplot.
#' 
#' @param slum A \code{slum} object. The experiment to plot the standard curves.
#' 
#' @seealso \code{\link{geom_sc}}
#' @author Renan Sauteraud
#' 
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point
#' @importFrom ggplot2 scale_x_log10 scale_y_log10
#' @export
plot_SC <- function(slum){
  msl <- melt(slum)
  msl_ss <- subset(msl, tolower(sample_type) == "standard")
  p <- ggplot(msl_ss, aes_string(color = "plate"), alpha = 0.5) + scale_x_log10() + 
    scale_y_log10() + facet_wrap(~analyte) + geom_sc(slum) + 
    geom_point(aes_string(x = "concentration", y = "mfi"))
  print(p)
}
