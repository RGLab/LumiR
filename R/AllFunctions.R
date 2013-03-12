#The root of the experiment
read.experiment<-function(path="./"){
  analyte.file<-list.files(path,pattern="analyte",full.names=TRUE)
  layout.file<-list.files(path,pattern="layout",full.names=TRUE)
  pheno.file<-list.files(path,pattern="phenotype",full.names=TRUE)
  ##TODO: throw error if some files are missing and tell the user to use setup_template
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
  phenoData<-.read.phenotype(path=path, pheno.file=pheno.file)

  if(length(layout.file)>0){
    layout<-.read.layout(layout.file)
    pData(phenoData)<-merge(pData(phenoData), layout, by="well")
  }

  #exprs
  if(type=="BIOPLEX"){
    phenoDT<-as.data.table(pData(phenoData))[,list(plate, well, sample_id)]
    exprs<-.read.exprs.bioplex(all.files)
    exprs<-merge(exprs, phenoDT, by=c("plate", "well"))
  } else {
    phenoDT<-as.data.table(pData(phenoData))[,list(plate, filename, sample_id)]
    if(type=="LXB"){ 
      exprs<-.read.exprs.lxb(all.files)
    } else {
      exprs<-.read.exprs.xPonent(all.files)
    }
    # Do we need to do this?
    setkeyv(phenoDT,c("plate","filename"))
    setkeyv(exprs,c("plate","filename"))    
    exprs<-exprs[phenoDT,]
  }
  setkey(exprs,sample_id)

  #fData
    featureData<-.read.analyte(analyte.file)
    mapping<-as.data.table(featureData@data)
    setkey(exprs,bid)
    setkey(mapping,bid)
    ## Join by bid
    exprs<-exprs[mapping,]
    # Re-order based on sample_id, then analyte
    setkeyv(exprs,c("sample_id","analyte"))
  
  blum<-new("blum", phenoData=phenoData, featureData=featureData, exprs=exprs)
  return(blum)
}

.getXponentWellsID<-function(filenames){
  filenames<-gsub(".csv", "", filenames)
  if(length(grep("Run", filenames))==0){#XPONENT v3.
    wellsID<-unlist(lapply(strsplit(filenames, split="_"), tail, 1))
  } else{#XPONENT v1.
    wellsNo<-as.numeric(gsub("Run", "", filenames))
    wellsID<-paste(LETTERS[(wellsNo-1)%%8+1],ceiling(wellsNo/8), sep="")
  }
  return(wellsID)
}
.getLXBWellsID<-function(filenames){#platename_wellID.lxb
  wellsID<-gsub(".lxb", "", unlist(lapply(strsplit(filenames, split="_"), tail, 1)))
  return(wellsID)
}
.getBioplexWellsID<-function(filenames){#file.xml > root>Wells
  wellsID<-unlist(lapply(filenames, function(x){
    root<-xmlRoot(xmlTreeParse(x))
    wNames<-xmlSApply(root[["Wells"]], function(x){
      paste(LETTERS[as.numeric(xmlAttrs(x)["RowNo"])], xmlAttrs(x)["ColNo"], sep="")
    })
    return(unlist(wNames))
  }))
  return(wellsID)
}

.read.exprs.xPonent<-function(filenames){
  exprsList<-lapply(filenames, function(x){dt<-fread(x);
                    ss=tail(strsplit(x,"/")[[1]],2)
                    fname=ss[2]  # Get plate & fname now for the merge with sample_ID
                    plate=ss[1]
                    dt[,c("plate","filename"):=list(plate,fname)]
                    #dt[,plate:=plate][,filename:=fname]
                    })
  ## rbind all data.tables
  ## The code could be improve when fread supports specifying the type for different columns
  exprs<-rbindlist(exprsList)
  exprs<-.sanitize.exprs(exprs)
  return(exprs)
}

.sanitize.exprs<-function(exprs){
  ## Lower case for all
  setnames(exprs,names(exprs),tolower(names(exprs)))
  bidGrep<-"id"
  bIdx<-grep(bidGrep, names(exprs))
  name<-names(exprs)  
  name[bIdx]<-"bid"
  setnames(exprs, name)
  ## Sanitize the names
  #namesToGrep<-"id|cl|dbl|dd|rp1|plate|name|time|well"  
  #cIdx<-grep(namesToGrep, names(exprs))
  #exprs<-exprs[,cIdx, with=FALSE]
  
  # Remove the eventno
  exprs<-exprs[,eventno:=NULL]  

  # Standardize the names
  setnames(exprs, "rp1", "fl")
  setnames(exprs, "cl1", "fl.x")
  setnames(exprs, "cl2", "fl.y")  
  
  ## Change the data types
  # intToGrep<-"id|cl|rp1|time"
  # iIdx<-grep(intToGrep, names(exprs))
  
  exprs<-exprs[,lapply(.SD, as.integer), by="plate,filename"]
#  exprs<-cbind(exprs[,names(exprs)[-iIdx], with=FALSE] ,exprs[,lapply(.SD, as.numeric), .SDcols=names(exprs)[iIdx]])
  return(exprs)
}

.read.exprs.lxb<-function(filenames){
  nFiles<-length(filenames)
  exprsList<-vector('list', nFiles)
  for(fileIdx in 1:nFiles){
    suppressWarnings(lxb<-read.FCS(filenames[[fileIdx]]))
    ss<-tail(strsplit(filenames[[fileIdx]],"/")[[1]],2)
    fname<-ss[2]  # Get plate & fname now for the merge with sample_ID
    plate<-ss[1]
    asdt<-as.data.table(exprs(lxb))
    asdt[,plate:=plate][,filename:=fname]
    exprsList[[fileIdx]]<-asdt
  }
  exprs<-rbindlist(exprsList)
  exprs<-.sanitize.exprs(exprs)
  return(exprs)
}
.read.exprs.bioplex<-function(filenames){
  exprsList<-list()
  for(filename in filenames){
    xml<-xmlTreeParse(filename)
    root<-xmlRoot(xml)
    exprsFile<-xmlApply(root[["Wells"]], function(x){
      plate<-tail(strsplit(filename,"/")[[1]],2)[1]
      well<-paste(LETTERS[as.numeric(xmlAttrs(x)["RowNo"])], xmlAttrs(x)["ColNo"], sep="")
      str<-xmlValue(x[["BeadEventData"]])
      ss<-unlist(strsplit(str, split="\n"))
      sss<-strsplit(ss, split=" ")
      dt<-as.data.table(do.call(rbind, lapply(sss, as.numeric)))
      setnames(dt, c("bid", "dd", "rp1", "cl1", "cl2"))
      dt<-dt[,plate:=plate][,well:=well]
      return(dt)
    })
    exprsList<-c(exprsList, exprsFile)
  }
  exprs<-rbindlist(exprsList)
  exprs<-.sanitize.exprs(exprs)
  return(exprs)
}

.read.analyte<-function(analyte.file){
  dt<-fread(analyte.file)
  setnames(dt, names(dt), tolower(names(dt)))
  if(length(dt)!=2 | !all(colnames(dt)%in%c("analyte","bid"))) {
    stop("The analyte mapping file should be a csv file with two columns 'analyte' and 'bid'\n")
  } else {
    featureData<-as(dt, "AnnotatedDataFrame")
  }
  return(featureData)
}

.read.layout<-function(layout.file){
  df<-read.csv(layout.file, header=TRUE)
  colnames(df)<-tolower(colnames(df))
  if(!all(colnames(df)%in%c("well", "sample_type", "concentration"))){#required cols
    stop("The layout mapping file should be a csv file with three columns 'well', 'sample_type' and 'concentration'\n")
  }
  if(length(unique(df$well))!=nrow(df)){#wells are unique
    stop("The layout file should contain only one line per well")
  }
  if(nrow(df[df$sample_type!="standard" & df$sample_type!="background" & !is.na(df$concentration),])){#conc only set for standards
    stop("The 'concentration' in layout mapping file should only be set for standard wells\n Check wells: ",paste(as.character(df[df$sample_type!="standard" & !is.na(df$concentration),"well"]), collapse=","))
  }
  if(nrow(df[df$sample_type=="background" & !(is.na(df$concentration) | df$concentration==0),])){
    stop("The 'concentration' for background samples in layout mapping file should be set to 0 or NA\n")
  }
  
  return(df)
}

.read.phenotype<-function(path, pheno.file){
  df <-read.csv(pheno.file, colClasses="factor")
  colnames(df)<-tolower(colnames(df))
  #df$concentration<-as.numeric(levels(df$concentration))[df$concentration] #More efficient than fact->char->numeric
  if(!all(c("plate","filename","well")%in%colnames(df))){
    stop("The phenotype mapping file must at least have the 'plate', 'filename' and 'well' columns\n")
  }
  for(i in 1:nrow(df)){
    if(length(list.files(paste(path, df[i,"plate"], sep="/"),pattern=as.character(df[i,"filename"])))==0){
      stop("The file ", as.character(df[i, "filename"]), " is not found in the given path. Verify plate and filename information in phenotype mapping file")
    }
  }
  phenoData<-as(df, "AnnotatedDataFrame")
  return(phenoData)
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
  featureData<-as(as.data.frame(mat), 'AnnotatedDataFrame')
  return(featureData)
}

### Summarize to MFIs and add standardCurves informations
slummarize<-function(from,type="MFI"){
  dt<-exprs(from)
  dt<-dt[,as.double(median(rp1)), by="sample_id,analyte"]
  setnames(dt, c("sample_id", "analyte", type))
  setkey(dt, sample_id)
  mat<-matrix(dt[,MFI], ncol=length(levels(dt[,sample_id])), dimnames=list(levels(dt[,analyte]),levels(dt[,sample_id])))
  #mat<-lapply(exprs(from),sapply,median)
  #mat<-t(do.call("rbind",mat))
  mfiSet<-new("slum", formula=as.formula("log(mfi) ~ c + (d - c)/(1 + exp(b * (log(x) - log(e))))^f"), inv=function(y, parmVec){exp(log(((parmVec[3] - parmVec[2])/(log(y) - parmVec[2]))^(1/parmVec[5]) - 1)/parmVec[1] + log(parmVec[4]))}
  )
  exprs(mfiSet)<-mat
  pData(mfiSet)<-pData(from)
  fData(mfiSet)<-fData(from)
  mfiSet@unit="MFI"
  
  df<-melt(mfiSet)
  # subselects standards
  df<-subset(df, concentration!=0 & tolower(sample_type)=="standard")		  
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
      conc<-c(conc, mfiSet@inv(mat[[i,j]], coefs[coefs$plate==strsplit(colnames(mat)[j], "_")[[1]][2] & coefs$analyte==rownames(mat)[i], 3:7]))
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

setup_templates<-function(path, templates=c("layout", "analyte", "phenotype")){
  
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
        fNames<-rep(unlist(lapply(all.files, lapply, function(x)tail(strsplit(x,"/")[[1]],1))), wellsPerFile)
        plate<-rep(plates, wellsPerFile)
        plate<-sapply(strsplit(plate, split="/"), tail, 1)
      } else {
        len<-lapply(plates, function(x){length(list_files_with_exts(x,exts=typeExt))})
        fNames<-unlist(lapply(all.files, function(x)tail(strsplit(x,"/")[[1]],1)))
        plate<-rep(plates, len)
        plate<-sapply(strsplit(plate, split="/"), tail, 1)
        if(type=="XPONENT"){
          wells<-.getXponentWellsID(fNames)
        } else if(type=="LXB"){
          wells<-.getLXBWellsID(fNames)
        } else {
          wells<-.getBioplexWellsID(fNames)
        }
      }
      noExt<-gsub(paste0(".",typeExt), "", fNames)
      sample_ID<-paste(noExt, plate, wells, sep="_")
      phenotype<-data.frame(plate=plate,filename=fNames,well=wells, sample_ID=sample_ID)
    
      write.csv(phenotype, file=paste0(path,"phenotype.csv"), row.names=FALSE)
    }
  }
    
  #layout
  if("layout"%in%templates){
    if(length(layout.file)>0){
      warning("The layout mapping file already exists, remove it to setup a template for it")
    } else {
      wells<-paste0(LETTERS[1:8], rep(seq(1,12), each=8))
      sample_type<-rep("unknown", length(wells))
      concentration<-rep(NA, length(wells))
      layout<-data.frame(well=wells, sample_type=sample_type, concentration=concentration)
      write.csv(layout, file=paste0(path,"layout.csv"), row.names=FALSE)
    }
  }

  #analyte
  if("analyte"%in%templates){
    if(length(analyte.file)>0){
      warning("The analyte mapping file already exists, remove it to setup a template for it")
    } else {
      if(type=="LXB" & length(list_files_with_exts(path, exts="lxd")>0)){
          lxdFile<-list_files_with_exts(path, exts="lxd")
          write.csv(pData(.read.lxd(lxdFile)), file=paste0(path, "analyte.csv"), row.names=FALSE)
      } else{
        if(type=="XPONENT"){  BIDs<-.getXponentBID(all.files[1])
        } else if(type=="LXB"){ BIDs<-.getLXBBID(all.files[1])
        } else { BIDs<-.getBioplexBID(all.files[1])
        }
        BIDs<-BIDs[BIDs!=0]
        analyte<-paste0(rep("unknown", length(BIDs)), BIDs)
        write.csv(data.frame(bid=BIDs, analyte=analyte), file=paste0(path,"analyte.csv"), row.names=FALSE)
      }
    }
  }
}

.getXponentBID<-function(firstFile){
  sLine<-grep("[Ee]vent[Nn]o", readLines(firstFile[1], n=5))-1
  con<-read.csv(firstFile, skip=sLine, header=TRUE);
  BIDs<-sort(as.numeric(unique(con[,2])))
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
  BIDs<-as.numeric(xmlSApply(root[["Wells"]][[1]][["RunSettings"]][["RegionsOfInterest"]], xmlAttrs))
  return(BIDs)
}
