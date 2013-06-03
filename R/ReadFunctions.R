#The root of the experiment
read.experiment<-function(path="./"){
  analyte.file<-list.files(path,pattern="analyte",full.names=TRUE)
  layout.file<-list.files(path,pattern="layout",full.names=TRUE)
  pheno.file<-list.files(path,pattern="phenotype",full.names=TRUE)
  missing_templates<-c()
  if(length(analyte.file)==0){
    missing_templates<-c(missing_templates, "analyte")
  }
  if(length(layout.file)==0){
    missing_templates<-c(missing_templates, "layout")
  }
  if(length(pheno.file)==0){
    missing_templates<-c(missing_templates, "phenotype")
  }
  if(!is.null(missing_templates)){
    cat(paste(missing_templates, collapse=", "), "mapping file missing.\nSome information may be missing.\n")
    template_list<-setup_templates(path, templates=missing_templates, write=FALSE)
  }

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
  if(length(pheno.file)>0){
    phenoData<-.read.phenotype(path=path, pheno.file=pheno.file)
  } else {
    phenoData<-as(template_list[["pheno"]], "AnnotatedDataFrame")
  }

  if(length(layout.file)>0){
    layout<-.read.layout(layout.file)
  } else {
    layout<-template_list[["layout"]]
  }
  pData(phenoData)<-merge(pData(phenoData), layout, by="well")

  #exprs
  phenoDT<-as.data.table(pData(phenoData))[,list(plate, filename, well, sample_id)]
  if(type=="BIOPLEX"){
    exprs<-.read.exprs.bioplex(all.files)
  } else if(type=="LXB"){ 
    exprs<-.read.exprs.lxb(all.files)
  } else {
    exprs<-.read.exprs.xPonent(all.files)
  }
  # Join with phenoData so we can add the sample_id to the exprs
  # Do we need to do this? Added well to make keys unique in phenoDT
  #setkeyv(phenoDT,c("plate","filename"))
  setkey(exprs,plate,filename,well)    
  exprs<-exprs[phenoDT,]
  setkey(exprs,sample_id)

  #fData
  if(length(analyte.file)>0){
    featureData<-.read.analyte(analyte.file)
  } else {
    featureData<-as(template_list[["analyte"]], "AnnotatedDataFrame")
  }
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
  filenames<-gsub(".csv", "", filenames, fixed=TRUE)
  filenames<-unlist(lapply(filenames, function(x)tail(strsplit(x,"/")[[1]],1)))
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

.sanitize.exprs<-function(exprs){
  setnames(exprs,names(exprs),tolower(names(exprs)))
  bidGrep<-"id"
  bIdx<-grep(bidGrep, names(exprs))
  name<-names(exprs)  
  name[bIdx]<-"bid"
  setnames(exprs, name)
  # Remove the eventno
  suppressWarnings(exprs<-exprs[,eventno:=NULL]) 
  # Standardize the names
  setnames(exprs, "rp1", "fl")
  setnames(exprs, "cl1", "fl.x")
  setnames(exprs, "cl2", "fl.y")  

  exprs<-exprs[,lapply(.SD, as.integer), by="plate,filename,well"]
  return(exprs)
}

.read.exprs.xPonent<-function(filenames){
  exprsList<-lapply(filenames, function(x){dt<-fread(x);
                    ss=tail(strsplit(x,"/")[[1]],2)
                    fname=ss[2]  # Get plate & fname now for the merge with sample_ID
                    plate=ss[1]
                    wname=.getXponentWellsID(fname)
                    dt[,c("plate","filename","well"):=list(plate,fname,wname)]
                    })
  ## rbind all data.tables
  ## The code could be improve when fread supports specifying the type for different columns
  exprs<-rbindlist(exprsList)
  exprs<-.sanitize.exprs(exprs)
  return(exprs)
}

.read.exprs.lxb<-function(filenames){
  nFiles<-length(filenames)
  wNames<-.getLXBWellsID(filenames)
  exprsList<-vector('list', nFiles)
  for(fileIdx in 1:nFiles){
    suppressWarnings(lxb<-read.FCS(filenames[[fileIdx]]))
    ss<-tail(strsplit(filenames[[fileIdx]],"/")[[1]],2)
    fname<-ss[2]  # Get plate & fname now for the merge with sample_ID
    plate<-ss[1]
    asdt<-as.data.table(exprs(lxb))
    asdt[,c("plate","filename","well"):=list(plate,fname,wNames[fileIdx])]
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
      fname<-tail(strsplit(filename,"/")[[1]],2)[2]
      well<-paste(LETTERS[as.numeric(xmlAttrs(x)["RowNo"])], xmlAttrs(x)["ColNo"], sep="")
      str<-xmlValue(x[["BeadEventData"]])
      ss<-unlist(strsplit(str, split="\n"))
      sss<-strsplit(ss, split=" ")
      dt<-as.data.table(do.call(rbind, lapply(sss, as.numeric)))
      setnames(dt, c("bid", "dd", "rp1", "cl1", "cl2"))
      dt<-dt[,c("plate","filename","well"):=list(plate,fname,well)]
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
    dt[,bid:=as.numeric(bid)]
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
    if(length(list.files(paste(path, df[i,"plate"], sep="/"),pattern=.strip.bad.filenames(as.character(df[i,"filename"]))))==0){
      stop("The file ", as.character(df[i, "filename"]), " is not found in the given path. Verify plate and filename information in phenotype mapping file")
    }
  }
  phenoData<-as(df, "AnnotatedDataFrame")
  return(phenoData)
}

.strip.bad.filenames<-function(string){
  #Escape the bad characters
  string<-gsub("\\+", "\\\\+", string)
  return(string)
}
