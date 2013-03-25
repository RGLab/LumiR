#
# saxcyb_utils function rewritten with less columns needed
#

#@mbama: A melted BAMAObject, or a BAMAObject?
#@trim: A logical, remove outliers if TRUE (5%)
splitAnalytes2<-function(mbama, trim=TRUE){
  #remove background and standards
  df<-subset(mbama, tolower(sample_type)!="standard" & tolower(sample_type)!="background" & tolower(sample_type)!="unmarked")
  df<-droplevels(df)
  analytes<-split(df, df$analyte)
  if(trim){
    analytes<-lapply(analytes, function(x){#one elt/bead
      trimidx<-unsplit(lapply(split(x, x$well), function(y){#one elt/well
        y$RP1 < mean(y$RP1,trim=0.05)+4*sd.trim(y$RP1,0.05)
        }), x$well)
        x[trimidx,]
      })
  }
  return(analytes)
}


sd.trim<-function(x, trim=0, na.rm=FALSE, ...){
# taken from http://ggorjan.blogspot.com/2008/11/trimmed-standard-deviation.html
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)){
    warning("argument is not numeric or logical: returning NA")
    return(NA_real_)
  }
  if(na.rm){x<-x[!is.na(x)]}
  if(!is.numeric(trim) || length(trim) != 1){
    stop("'trim' must be numeric of length one")}
  n<-length(x)
  if(trim>0 && n>0){
    if(is.complex(x)) {stop("trimmed sd are not defined for complex data")}
    if(trim >= 0.5) {return(0)}
    lo<-floor(n*trim)+1
    hi<-n+1-lo
    x<-sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  return(sd(x))
}

#INPUT: data.frame
#OUTPUT: data.frame with added column for repeat_effect
#@analyte: data.frame containing every bead for the analyte
#@mbpg: min bead per group (all repeats added)
#@mbpw: min bead per well
adjustRepeats2<-function(analyte, mbpg=99, mbpw=30){
# modify the analyte data frame so that wells with insufficient number of beads are merged
# estimate repeats effects
  an_group<-split(analyte, analyte$sample_id)
  pvals <- numeric()
  ret <- data.frame()

  for ( grp in 1:length(an_group) ) {
    if(nrow(an_group[[grp]]) > mbpg) {
      an_group[[grp]]$well <- drop.levels(an_group[[grp]]$well)
      numbeads<-summary(an_group[[grp]]$well)    #number of beads per rep
      if(min(numbeads) > mbpw && length(levels(an_group[[grp]]$well)) > 1) {
        #Huber regression
        fit_huber <- rlm(RP1~well,data=an_group[[grp]],contrasts=list(well="contr.sum"),scale.est='proposal 2',k=1.345,maxit=100)
        coef.rep <- c(coef(fit_huber)[-1],-sum(coef(fit_huber)[-1]))
        names(coef.rep) <- levels(an_group[[grp]]$well)
        for ( k in levels(an_group[[grp]]$well) ) { 
          an_group[[grp]][an_group[[grp]]$well==k,"rep_effect"] <- coef.rep[k] #adds a column for rep_effect
        }   
      } else {     # multiple reps but insufficient beads
        cat('too few beads per well: group_name=', names(an_group)[grp], ', beads=', nrow(an_group[[grp]]), ', well=', which.min(numbeads), '(',min(numbeads),'). merged with other wells\n',sep="")
        an_group[[grp]][,"rep_effect"]<-0
        #for ( k in levels(an_group[[grp]]$well) ) { 
          #an_group[[grp]][an_group[[grp]]$well==k,"rep_effect"] <- 0
        #}   
      }   
      ret<-rbind(ret, an_group[[grp]])
    } else { # too few pooled beads
      cat('too few pooled beads: group_idx=', names(an_group)[grp], ', beads=', nrow(an_group[[grp]]), '. skip.\n',sep="")
    }   
  }   
  return(ret)
}

logtransfit.group2 <- function( analyte, SB=1, maxIter=5 ) {
# logtrans (library MASS) regression fitting for SAxCyB ANOVA model
# - simultaneously fit independent groups
  ctrl_lvls <- levels(analyte$control_idx)
  # construct a regession formula
  contrlist <- list(group_name="contr.sum")
  f <- as.formula("y~group_name")
  y <- analyte$RP1 - analyte$rep_effect - SB
  shifts <- logtrans(f, data=analyte, contrasts=contrlist, plotit=FALSE )
  best_shift <- shifts$x[which.max(shifts$y)[1]]
  y <- log( y + best_shift )
  analyte$y <- y

  ret<-vector('list', length(ctrl_lvls))
  names(ret)<-ctrl_lvls
  for(ctrllvl in ctrl_lvls){
    cur_ctrlgrp<-subset(analyte, control_idx==ctrllvl) #a big d.f with all the groups matching the control
    cur_ctrlgrp<-droplevels(cur_ctrlgrp)
    ctrl_name<-as.character(unique(cur_ctrlgrp[cur_ctrlgrp$sample_type=="control", "group_name"]))
    cur_ctrlgrp$group_name<-relevel(cur_ctrlgrp$group_name, ref=ctrl_name) #control first
    if(length(levels(cur_ctrlgrp$group_name))<2) {next}
    for(iter in 1:maxIter){
      cur_fit0<-lm(f, data=cur_ctrlgrp)
      w<-1/unsplit(lapply(split(resid(cur_fit0),cur_ctrlgrp$group_name), function(x){rep(mean(x^2),length(x))}),  cur_ctrlgrp$group_name)
      cur_fit <- lm(f, data=cur_ctrlgrp, weights=w )
      upd<-updateRepeatEffects2(cur_ctrlgrp, cur_fit)
      cur_ctrlgrp$rep_effect<-cur_ctrlgrp$rep_effect+upd$expanded_delta_rep_effect
      D<-unlist(upd$delta_rep_effect)
      if(mean(D^2)<0.1) {break}
    }
    nGrps<-length(levels(cur_ctrlgrp$group_name))-1 #number of groups associated with the current ctrl
    h_names<-as.character(unique(cur_ctrlgrp[cur_ctrlgrp$sample_type!="control", "group_name"]))
    K<-cbind(rep(0,nGrps), diag(1, nGrps, nGrps))
    rownames(K)<-h_names
    #General linear hypothesis multiple comparisons
    cur_ht<-glht(cur_fit, linfct=K)
    cur_coef<-coef(cur_fit)
    ret[[ctrllvl]]<-list(ht=cur_ht,coef_full=cur_coef,shift=best_shift,ctrllvl=ctrllvl,contrasts='treatment')
  }
  rets=list(fits=ret, transformedFI=y)
  return(rets)
}



#@cur_ctrlgrp: A data.frame with all the groups of a given control_idx
#@cur_fit: A lm object
#OUTPUT: A vector of numeric to update the rep_effect (beta in the paper) col of the d.f
updateRepeatEffects2<-function(cur_ctrlgrp, cur_fit){
  an_group<-split(cur_ctrlgrp, cur_ctrlgrp$group_name)
  gam<-coef(cur_fit)[1]+c(0, coef(cur_fit)[-1]) #fitted means
  delta_rep_effect<-vector('list', length(an_group))
  expanded_delta_rep_effect<-numeric()
  for(i in 1:length(an_group)){
    curGrp<-an_group[[i]]
    curGrp<-droplevels(curGrp)
    #gradient
    g<-unlist(lapply(split(curGrp, curGrp$well),function(x){
      -sum((x$y-gam[i])/exp(x$y))}))
    #Hessian
    H<-unlist(lapply(split(curGrp, curGrp$well),function(x){
      sum((x$y+1)/exp(2*x$y))}))
    if(all(H>1e-5)){ #equality constrained Newton
      descent_dir<- -g/H - sum(-g/H)/H/(sum(1/H))
    } else{ #projected gradient
      descent_dir<- -g+mean(g)
    }
    #expand descent direction for evaluating the objective
    names(descent_dir)<-levels(curGrp$well)
    exp_descent_dir<-numeric(nrow(curGrp))
    for(kk in levels(curGrp$well)){
      exp_descent_dir[curGrp$well==kk]<-descent_dir[kk]
    }
    #backtracking line search
    t<-1
    alpha_param<-0.25
    beta_param<-0.8
    cur_objective<-sum((curGrp$y-gam[i])^2)
    for(iter in 1:10){
      new_term<-exp(curGrp$y)-t*exp_descent_dir
      if(all(new_term>1e-5)){
        new_objective<-sum((log(new_term)-gam[i])^2)
        if(new_objective-alpha_param*t*sum(g*descent_dir)<cur_objective) {break}
      }
      t<-beta_param*t
    }
    delta_rep_effect[[i]]<-t*descent_dir
    expanded_delta_rep_effect<-c(expanded_delta_rep_effect,t*exp_descent_dir)
  }
  ret<-list(delta_rep_effect=delta_rep_effect, expanded_delta_rep_effect=expanded_delta_rep_effect)
  return(ret)
}
    
#@analyte: A data.frame containing a single analyte
#@mSB: A numeric, the mean intensity of the trimmed beads for this analyte
saxcyBfit<-function(analyte, mSB){
  an_repEffect<-adjustRepeats2(analyte)
  model<-logtransfit.group2(an_repEffect, mSB)
  fits<-model$fits
  an_repEffect$transformedFI<-model$transformedFI
  return(list(fits=fits, analyte=an_repEffect))
}


#@analyte: A data.frame containing a single analyte
#OUTPUT: A data.frame with the MFI in each well for the analyte
getMFI2<-function(analyte){
  ret<-aggregate(analyte$RP1, list(
    group_name=analyte$group_name,
    sample_type=analyte$sample_type,
    control_idx=analyte$control_idx,
    well=analyte$well), median)
  names(ret)[5]="MFI"
  return(ret)
}

SAxCyBpower.ineq <- function( sfit, Th, effectsize, alpha=0.05 ) { 
# estimate power for SAxCyB bioinequivalence test
#  input        sfit= a SAxCyB fit
#               Th= equivalence margin
#               effectsize= effect size at which to estimate the power
#               alpha= desired significance level
#  output       estimated power
  ht <- sfit$ht
  z_alpha <- qnorm(1-alpha)
  delta_L = (effectsize+Th)/sqrt(diag(vcov(ht)))
  delta_R = (effectsize-Th)/sqrt(diag(vcov(ht)))

  power.ineq=array(length(delta_L))
  for ( i in 1:length(delta_L) ) { 
    power.ineq[i]=1-pmvnorm(lower=c(z_alpha,-Inf),upper=c(Inf,-z_alpha),mean=c(delta_L[i],delta_R[i]))
  }
  return(power.ineq)
}

SAxCyBpval <- function( sfit, Th ) {
# p-value computation for SAxCyB model
#  input        sfit= a SAxCyB fit
#               Th= equivalence margin 
#  output       p.schuirmann_t= vector of p-values of a set of Schuirmann-type bioinequivalence tests ( two one-sided t-tests)
#               p.schuirmann_t.mcp= p.schuirmann_t adjusted for multiple hypothesis testing
  ht <- sfit$ht
  zval<-coef(ht)/sqrt(diag(vcov(ht)))
  # univariate
  t_over <- zval + Th/sqrt(diag(vcov(ht)))
  t_under <- zval - Th/sqrt(diag(vcov(ht)))
  
  p1 <- pnorm(t_over)
  p2 <- 1-pnorm(t_under)
  p.schuirmann_t <- pmin(p1,p2)
  p.schuirmann_t.mcp <- p.adjust(p.schuirmann_t,method="holm")
  return(list( p.schuirmann_t=p.schuirmann_t, p.schuirmann_t.mcp=p.schuirmann_t.mcp ) ) 
}

simplepval <- function( fit, MFI, analyte ) { 
## wrapper function for MFIttestpval2 and fullFIttestpval2
# simple MFI-based t-test
  p.MFI_t <- MFIttestpval2( MFI ) 
  p.MFI_t <- p.MFI_t[rownames(fit$ht$linfct)]
  p.MFI_t.mcp <- p.adjust(p.MFI_t, method="holm")
  
  # full FI-based t-test
  p.fullFI_t <- fullFIttestpval2( analyte ) 
  p.fullFI_t <- p.fullFI_t[rownames(fit$ht$linfct)]
  p.fullFI_t.mcp <- p.adjust(p.fullFI_t, method="holm")
  
  return(list( p.MFI_t=p.MFI_t, p.MFI_t.mcp=p.MFI_t.mcp, p.fullFI_t=p.fullFI_t, p.fullFI_t.mcp=p.fullFI_t.mcp ) )
}

MFIttestpval2 <- function( MFI ) {
# t-test based on median FI
  cMFI<-MFI[MFI$sample_type=="control","MFI"]
  cases<-MFI[MFI$sample_type!="control",]
  cases<-droplevels(cases)
  z<-split(cases, cases$group_name)
  samples<-names(z) #group_names
  sName <- character()
  p.ttest <- numeric()
  for ( sample in samples ) {
    sMFI <- z[[sample]]$MFI
    p.ttest <- c( p.ttest, tryCatch(t.test(cMFI,sMFI)$p.value,error=function(e) {NA} ) )
    sName <- c( sName, as.character(z[[sample]]$group_name[1]) )
  }
  names(p.ttest) <- sName
  return(p.ttest)
}

fullFIttestpval2 <- function( analyte ) {
# t-test based on all FI measurements
  cFI<-analyte[analyte$sample_type=="control","RP1"]
  cases<-analyte[analyte$sample_type!="control",]
  cases<-droplevels(cases)
  z<-split(cases, cases$group_name)
  samples<-names(z)
  sName <- character()
  p.ttest <- numeric()
  for ( sample in samples ) {
    sFI <- z[[sample]]$RP1
    p.ttest <- c( p.ttest, tryCatch(t.test(cFI,sFI)$p.value,error=function(e) {NA} ) )
    sName <- c( sName, as.character(z[[sample]]$group_name[1]) )
  }
  names(p.ttest) <- sName
  return(p.ttest)
}

fullFIttestpval <- function( analyte ) {
# t-test based on all FI measurements
  z<-split(analyte,analyte$group_name)
  sample_indices <- names(z)
  cFI <- z[["control"]]$RP1       # pool repeat wells
  sName <- character()
  p.ttest <- numeric()
  for ( group_name in sample_indices ) {
    if ( group_name != 'control' && nrow(z[[group_name]])>0 ) {
      sFI <- z[[group_name]]$RP1
      p.ttest <- c( p.ttest, tryCatch(t.test(cFI,sFI)$p.value,error=function(e) {NA} ) )
      sName <- c( sName, as.character(z[[group_name]]$group_name[1]) )
    }
  }
  names(p.ttest) <- sName
  return(p.ttest)
}

selThreshold<-function(deltaList, powerThresholds=c(0.19,0.29,0.39,0.49,0.59,0.69,0.79)){
  selDelta<-sapply(deltaList, median, na.rm=TRUE)
  curvature <- diff(diff(selDelta))
  signs <- sign(curvature)
  first_sign <- signs[1]
  inflection_idx <- min(which(signs!=first_sign))
  if ( is.finite(inflection_idx) ) {
    threshold <- powerThresholds[inflection_idx+1]
  } else { # if no inflection point found, find the point of the largest drop
    threshold <- powerThresholds[which.min(diff(selDelta))+1]
  }
  return(threshold)
}

#R#Probably no need for this function / or maybe just once for the actual PowerThreshold
popCSV2<-function(filename, isAppended, fit, MFI, analyte, pvals, delta, mSB, sdSB){
  bead_id = unique(analyte$bid)
  bead_name = unique(analyte$analyte)
  #^done
  rep_eff<-unlist(lapply(split(analyte,analyte$well),function(x) {x$rep_effect[1]}))
  b<-data.frame(well=levels(analyte$well),beta=rep_eff[levels(analyte$well)])
  mb<-merge(mfi,b)
  betastr <- unlist(lapply(split(mb,mb$group_name), function(x) { paste(format(x$beta,digits=5,trim=T),collapse=";")} ))
  num_reps0 <- unlist(lapply((tapply(analyte$rep, analyte$group_idx, unique)),length))
  grp_names <- levels(analyte$group_idx)
  exp_names <- as.character(unlist(lapply(split(analyte,analyte$group_idx),function(x) {x$group_name[1]} )))
  well <- unlist(lapply(split(MFI,MFI$group_name), function(x) { paste(x$well,collapse=";")} ))
  MFIstr <- unlist(lapply(split(MFI,MFI$group_name), function(x) { paste(x$MFI,collapse=";")} ))
  sample_idx <- unlist(lapply(split(analyte,analyte$group_name), function(x) { x$sample_idx[1] } ))
  control_idx <- unlist(lapply(split(analyte,analyte$group_name), function(x) { x$control_idx[1] } ))
  group_idx <- unlist(lapply(split(analyte,analyte$group_name), function(x) { x$group_idx[1] } ))

        if ( is.null(fit$contrasts) ) { # sum constraint
                alphas <- fit$coef[2:length(grp_names)]
                alphas <- c( alphas, -sum(alphas) )
        } else {        # treatment contrast
                alphas <- fit$coef[2:length(grp_names)]
                alphas <- c( 0, alphas )
        }
        names(alphas) <- names(sort(group_idx))
        diff.subj <- coef(fit$ht)

        sd.diff.subj <-  sqrt(diag(vcov(fit$ht)))
        p.diff.uni <- pvals$p.schuirmann_t
        sigcode.diff.uni <- unlist(lapply(p.diff.uni, sigcode))
        p.diff.uni.mcp <-  pvals$p.schuirmann_t.mcp
        sigcode.diff.uni.mcp <- unlist(lapply(p.diff.uni.mcp, sigcode))
        # t-test based on MFIs
        p.MFIttest <- pvals$p.MFI_t
        p.MFIttest<-p.MFIttest[names(p.diff.uni)]
        sigcode.MFIttest <- unlist(lapply(p.MFIttest, sigcode))
        p.MFIttest.mcp <- pvals$p.MFI_t.mcp
        sigcode.MFIttest.mcp <- unlist(lapply(p.MFIttest.mcp, sigcode))
        # t-test based on full FIs
        p.fullFIttest <- pvals$p.fullFI_t
        p.fullFIttest<-p.fullFIttest[names(p.diff.uni)]
        sigcode.fullFIttest <- unlist(lapply(p.fullFIttest, sigcode))
        p.fullFIttest.mcp <- pvals$p.fullFI_t.mcp
        sigcode.fullFIttest.mcp <- unlist(lapply(p.fullFIttest.mcp, sigcode))

        maxMFI <- aggregate(MFI$MFI,list(group_name=MFI$group_name), max)
        names(maxMFI)[2] <- "maxMFI"
        nreps <- num_reps0
        names(nreps) <- exp_names
        nreps<- nreps[as.character(maxMFI$group_name)]
        basic <- data.frame( bead_id=bead_id, bead_name=bead_name, group_name=maxMFI$group_name, sample_idx=sample_idx[as.character(maxMFI$group_name)], control_idx=control_idx[as.character(maxMFI$group_name)], well=well, rep=nreps,  MFI=MFIstr, maxMFI=maxMFI$maxMFI, mSB=mSB, sdSB=sdSB, delta=delta, alpha=alphas[as.character(maxMFI$group_name)], beta=betastr )

        ma<-data.frame( bead_id=bead_id, bead_name=bead_name, sample_idx=sample_idx[as.character(names(diff.subj))], control_idx=control_idx[as.character(names(diff.subj))], group_name=names(diff.subj), diff.subj=diff.subj, sd.diff.subj=sd.diff.subj, p.diff.uni=p.diff.uni, sig.p.diff.uni=sigcode.diff.uni, p.diff.uni.mcp=p.diff.uni.mcp, sig.p.diff.uni.mcp=sigcode.diff.uni.mcp, p.MFIttest=p.MFIttest, sig.p.MFIttest=sigcode.MFIttest, p.MFIttest.mcp=p.MFIttest.mcp, sig.p.MFIttest.mcp=sigcode.MFIttest.mcp, p.fullFIttest=p.fullFIttest, sig.p.fullFIttest=sigcode.fullFIttest, p.fullFIttest.mcp=p.fullFIttest.mcp, sig.p.fullFIttest.mcp=sigcode.fullFIttest.mcp )

        mtot <- merge(basic,ma,by=c("bead_id","bead_name","control_idx","sample_idx","group_name"),all=T)

        write.table( mtot, file=filename, append=isAppended, row.names=FALSE, col.names=!isAppended, sep=",", quote=FALSE )
}

saxcyB<-function(bama){
	mbama<-melt(bama)
	san<-splitAnalytes2(mbama)
	mSB <- lapply(san, function(x){mean(x$RP1, trim=0.05)})
	sdSB<- lapply(san, function(x){sd.trim(x$RP1, trim=0.05)})
	deltas <- c(0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32)
	PowerThresholds <- c(0.19,0.29,0.39,0.49,0.59,0.69,0.79)
	
	pList<-vector('list', length(PowerThresholds))
	alist<-beadsList<-ctrList<-groups<-c()
	
	for(idxk in 1:length(san)){
		saxcyObj <- saxcyBfit( san[[idxk]], mSB[[idxk]] )
		mfi2 <- getMFI2( saxcyObj$analyte )
		
		
		for ( idxFit in 1:length(saxcyObj$fits) ) {
			# select delta by power estimation
			power.est<-array(length(deltas))
			for ( delta.idx in 1:length(deltas) ) {
				power.est[delta.idx]<-mean(SAxCyBpower.ineq(saxcyObj$fits[[idxFit]], deltas[delta.idx], abs(coef(saxcyObj$fits[[idxFit]]$ht)) ))
			}
			mfi <- drop.levels(subset(mfi2,control_idx==saxcyObj$fits[[idxFit]]$ctrllvl) )
			an <- subset(saxcyObj$analyte,control_idx==saxcyObj$fits[[idxFit]]$ctrllvl )
			
			# compact group ids
			an$group_name = drop.levels(an$group_name)
			ctrl_name <- as.character(unique(an[an$sample_type=="control", "group_name"]))
			
			alpha<-saxcyObj$fits[[idxFit]]$coef_full
			len<-length(alpha)-1 #groups in the fit - ctrl
			groups<-c(groups, substr(names(alpha[-1]), 11, nchar(names(alpha[-1]))))
			ctrList<-c(ctrList, rep(ctrl_name, len))
			aList<-c(aList,alpha[-1])
			beadsList<-c(beadsList, rep(names(san)[[idxk]], len))
			
			# reset reference group id to control
			an$group_name<-relevel(drop.levels(an$group_name),ref=ctrl_name)			
			qvals<-simplepval(saxcyObj$fits[[idxFit]], mfi, an )
			
			#save the delta and pvals
			for ( PowerThreshold in PowerThresholds ) {
				pt<-as.character(PowerThreshold)
				delta<-deltas[max(which(power.est>PowerThreshold))]
				dList[[pt]]<-c(dList[[pt]], rep(delta, length(levels(an$group_name))))
				pList[[pt]]<-c(pList[[pt]], SAxCyBpval(saxcyObj$fits[[idxFit]], delta)$p.schuirmann_t.mcp)
				
			}
		}
	}
	thresh<-selThreshold(dList, PowerThresholds) #Select threshold using deltas
	df<-data.frame(bead=beadsList, group=groups, control=ctrList, alpha=aList, p.val=pList[[as.character(thresh)]])
	return(df)
}
