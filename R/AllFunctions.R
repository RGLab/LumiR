####################
##    TEMPLATES   ##
####################
setup_templates<-function(path, templates=c("layout", "analyte", "phenotype"), write=TRUE){
  dfList<-list()
  analyte.file<-list.files(path,pattern="analyte",full.names=TRUE)
  layout.file<-list.files(path,pattern="layout",full.names=TRUE)
  pheno.file<-list.files(path,pattern="phenotype",full.names=TRUE)
  plates<-list.dirs(path, recursive=FALSE)

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
        #plate<-sapply(strsplit(plate, split="/"), tail, 1)
        plate <- basename(plate)
      } else {
        len<-lapply(plates, function(x){length(list_files_with_exts(x,exts=typeExt))})
        #fNames<-unlist(lapply(all.files, function(x)tail(strsplit(x,"/")[[1]],1)))
        fNames <- basename(all.files)
        plate<-rep(plates, len)
        #plate<-sapply(strsplit(plate, split="/"), tail, 1)
        plate <- basename(plate)
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
      phenotype.df<-data.frame(plate=plate,filename=fNames,well=wells, sample_id=sample_id)
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
          analyte.df<-.read.lxd(lxdFile) 
      } else if(type=="BIOPLEX"){
        analyte.df<-.getBioplexAnalytes(all.files[1])
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
  sample_type <- tolower(unlist(xmlSApply(node, names)))
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

#Make results similar to https://immport.niaid.nih.gov/example_submission_packages/MBAA_Results.xls
#Source ID, Source ID Type,  Assay ID,  Assay Group ID,  Analyte Name,  MFI, Concentration Value, Concentration Unit,  MFI, Coordinate,  Comment
writeMBAA <- function(object, outfile="./MBAA_results.csv", concentration_unit="pg/mL"){
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
  write.csv(OUT, file=outfile, row.names=FALSE)
  #return(OUT)
}


####################
##  SUMMARY       ##
####################
### Summarize to MFIs and add standardCurves informations
slummarize<-function(from,type="MFI"){
  dt<-exprs(from)
  dt<-dt[,as.double(median(fl)), by="sample_id,analyte"]
  setnames(dt, c("sample_id", "analyte", type))
  setkey(dt, sample_id)
  mat<-matrix(dt[,MFI], ncol=length(unique(dt[,sample_id])), dimnames=list(unique(dt[,analyte]),unique(dt[,sample_id])))
  mfiSet<-new("slum", formula=as.formula("log(mfi) ~ c + (d - c)/(1 + exp(b * (log(x) - log(e))))^f"), inv=function(y, parmVec){exp(log(((parmVec[3] - parmVec[2])/(log(y) - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + log(parmVec[4]))}
  )
  exprs(mfiSet)<-mat
  pData(mfiSet)<-pData(from)
  fData(mfiSet)<-fData(from)
  mfiSet@unit="MFI"

  df<-melt(mfiSet)
  # subselects standards
  df<-subset(df, concentration!=0 & tolower(sample_type)=="standard")
  if(nrow(df)==0){
    stop("The object does not contain any standard. Edit layout.csv to include standards location and concentration.")
  }
  # Split by plate
  sdf<-split(df,df$plate)
  df2<-lapply(sdf,.fit_sc, mfiSet@inv)
  df2<-do.call("rbind",df2)
  mfiSet@fit<-df2

  conc<-c()
  coefs<-unique(mfiSet@fit[,c("plate", "analyte", "b","c","d","e","f")])
  inv<-mfiSet@inv
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      conc<-c(conc, as.numeric(mfiSet@inv(mat[[i,j]],
               coefs[coefs$plate==pData(mfiSet)[pData(mfiSet)$sample_id==colnames(mat)[j], "plate"] & coefs$analyte==rownames(mat)[i], 3:7])))
    }
  }
  concMat<-matrix(conc, ncol=ncol(mat))
  rownames(concMat)<-rownames(mat)
  colnames(concMat)<-colnames(mat)
  assayData(mfiSet)<-list(exprs=mat, concentration=concMat)
  mfiSet
}

.fit_sc<-function(df, inv)
{
  nCtrl<-length(unique(df$well)) #number of wells with standards
  df.split<-split(df, df$analyte)
  coeffs<-lapply(df.split, function(x){
    res<-drm(log(mfi) ~ concentration, data=x,fct=LL.5())
    return(res$parmMat)
  })

  calc_conc<-p100rec<-numeric(nrow(df))
  for(idx in 1:nrow(df))
  {
    calc_conc[idx]<-inv(df[idx,"mfi"], coeffs[[df[idx,"analyte"]]])
    p100rec[idx]<-calc_conc[idx]/df[idx,"concentration"]*100
  }
  sortCoeffs<-do.call("rbind", lapply(coeffs[df$analyte], t))
  colnames(sortCoeffs)<-c('b','c','d','e','f')
  df2<-cbind(df[,c("sample_id", "plate", "filename", "well", "analyte", "mfi", "concentration")], calc_conc, p100rec, sortCoeffs)
  return(df2)
}
