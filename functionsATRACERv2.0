##################################################################################
##################################################################################

# UTILITY FUNCTIONS FOR ALCA

##################################################################################
##################################################################################


#  factorList
#=================================================================================

# Determine the number of factors and treatments
factorlist<-function(inputparameters=inputparameters){
  #Define the input variables
  folder.output=inputparameters$folder.output
  condlist<-inputparameters$condList
  nfactor<-inputparameters$nFactor
  factorNames<-inputparameters$namesFactors
  nexp<-inputparameters$nexp
  totalfiles<-length(inputparameters$filesFluc)
  tfloc<-inputparameters$tfloc
  nplate<-inputparameters$nplate
  if (tfloc=="RAND"){
  	nfpexp<-totalfiles/(nexp*nplate)
  } else {
	nfpexp<-totalfiles/(nexp)
  }
  timevector<-inputparameters$timeVector
  ntime=length(timevector)
  inputparameters$ntime<-ntime
  inputparameters$nfpexp<-nfpexp

  expList<-NULL
  for (aa in 1:length(condlist)){
	expList<-rbind(expList,data.frame(read.delim(condlist[aa], header = FALSE)))
  }
  if (nfactor==1){
	TFIds <- levels(expList[,2]) 
	TFIds <-TFIds[order(TFIds)]
	treatmentIds <- levels(expList[,3])      
	treatmentIds <-treatmentIds[order(treatmentIds)]
	factorIds<-namesFactors
	levelsIds<-list()
	levelsIds[[factorIds]]<-treatmentIds
	expList<-data.frame(expList,expList[,3])
      colnames(expList)<-c("ID","TF","Treatment",factorIds)
  } else {
	expList2<-expList
	expList<-NULL
	expListtemp2<-apply(expList2, 1, function(x){
	  temp<-as.character(x[3])
	  for (i in 4:length(x)){
	  	temp<-paste(temp,as.character(x[i]),sep=".")
	  }
	  return(temp)
	})
	TFIds <- levels(expList2[,2])
	TFIds <-TFIds[order(TFIds)]
	treatmentIds <- levels(as.factor(expListtemp2))
	treatmentIds <-treatmentIds[order(treatmentIds)]
	levelsIds<-as.list(apply(expList2[,3:ncol(expList2)],2,function(x){levels(as.factor(x))}))
	for (i in 1:length(factorNames)){
		levelsIds[[factorNames[i]]]<-levelsIds[[i]][order(levelsIds[[i]])]
	}
	factorIds<-factorNames
	expList<-data.frame(expList2[,1:2],as.factor(expListtemp2),expList2[,3:ncol(expList2)])
      colnames(expList)<-c("ID","TF","Treatment",factorIds)
  }
  ExpCond<-list(tfNames=TFIds,treatmentNames=treatmentIds,levelsNames=levelsIds,
	factorNames=factorIds,expList=expList)
  ntreat<-length(treatmentIds)
  inputparameters$ntreat<-ntreat
  exptemp<-list(inputparameters=inputparameters, ExpCond=ExpCond)
  return(exptemp)
}

# readData
#=============================================================================================

readData<-function(dataIn=dataIn,mapIn=mapIn,expListIn=expListIn,plateIn=nplate,
  	inputparameters=inputparameters){
	  folder.output=inputparameters$folder.output
  if (nplate>1){
	datatempmp<-NULL
	maptempmp<-NULL
      listtempmp<-NULL
	expListtempm<-NULL
	for (gg in 1:plateIn){
		datatemp<-read.delim(dataIn[gg], header = FALSE)
		pos_first<-min(which(apply(datatemp,c(1,2),function(x) !is.na(as.numeric(x)))=="TRUE"))# The warning message for NA is ok
		col_first<-trunc(pos_first/dim(datatemp)[1])+1
		row_first<-pos_first-(col_first-1)*dim(datatemp)[1]
		datatemp <- datatemp[row_first:dim(datatemp)[1], col_first:dim(datatemp)[2]]
		datatemp <- as.numeric(as.matrix(datatemp))
		datatempmp<-c(datatempmp,datatemp)
		maptemp<-read.delim(mapIn[gg], header = FALSE)
		maptemp<-as.numeric(as.matrix(maptemp))
		maptempmp<-c(maptempmp,maptemp)
		listtemp<-read.delim(expListIn[gg], header = FALSE)
            listtempmp<-rbind(listtempmp,listtemp)
	}
	maptemp<-maptempmp
	datatemp<-datatempmp
	expListtemp<-listtempmp
  }else{
	datatemp<-read.delim(dataIn, header = FALSE)
	pos_first<-min(which(apply(datatemp,c(1,2),function(x) !is.na(as.numeric(x)))=="TRUE"))# The warning message for NA is ok
	col_first<-trunc(pos_first/dim(datatemp)[1])+1
	row_first<-pos_first-(col_first-1)*dim(datatemp)[1]
	datatemp <- datatemp[row_first:dim(datatemp)[1], col_first:dim(datatemp)[2]]
	datatemp <- as.numeric(as.matrix(datatemp))
	maptemp<-as.numeric(as.matrix(read.delim(mapIn, header = FALSE)))
	expListtemp<-read.delim(expListIn, header = FALSE)
  }
  dataInput<-list(data=datatemp,map=maptemp,expList=expListtemp)
  return(dataInput)
}

# fillSummary
#=========================================================================================

fillSummary<-function(counter=files,inputparameters=inputparameters,datatable=dataAll,
  expCond=expCond,expListtemp=expListIn){
    folder.output=inputparameters$folder.output
  treatIdstemp<-datatable$treatments
  ntreattemp<-length(treatIdstemp)
  ntreat<-inputparameters$ntreat
  nplate<-inputparameters$nplate
  treatIds<-expCond$treatmentNames
  timevector<-inputparameters$timeVector
  ntime<-inputparameters$ntime
  tfids<-expCond$tfNames
  ntf<-length(tfids)
  nfpexp<-inputparameters$nfpexp
  tfidstemp<-unique(datatable$expList[,2])
  if (nplate==1){
  	summary<-data.frame(matrix(NA,nrow=ntreattemp,ncol=(ntf+6)))
  	summary[,1]<-ceiling(counter/nfpexp)
  	summary[,2]<-ifelse(counter%%nplate==0,nplate-counter%%nplate,
		counter%%nplate)
  	summary[,3]<-treatIdstemp 
  	summary[,4]<-timevector[ifelse(counter%%(nplate*ntime)==0,
		(ntime*nplate-counter%%(ntime*nplate)),counter%%(ntime*nplate))]
  	for (i in 1:ntf){
  		summary[,4+i]<-ifelse(is.na(match(tfids,tfidstemp)[i]),0,1)
  	}
  } else {
	summary<-data.frame(matrix(NA,nrow=(ntreattemp*nplate),ncol=(ntf+6)))
	summary[,1]<-ceiling(counter/nfpexp)
  	summary[,2]<-rep(seq(1,nplate,by=1),each=ntreattemp)
	summary[,3]<-treatIdstemp 
  	summary[,4]<-timevector[ifelse(counter%%ntime==0,
		(ntime-counter%%ntime),counter%%ntime)]
      for (hh in 1:nplate){
  		listtemp<-read.delim(expListtemp[hh], header = FALSE)
		tfidstemp<-unique(listtemp[,2])
  		for (i in 1:ntf){
  			summary[(1+(hh-1)*ntreattemp):(ntreattemp*hh),4+i]<-
				ifelse(is.na(match(tfids,tfidstemp)[i]),0,1)
		}
  	}
  }	
  return(summary)
}

# orderTables
#=====================================================================================
 orderTables<-function(datatable=dataAll,expCond=expCond,inputparameters=inputparameters,
 	counter=files,index=bb){
	  folder.output=inputparameters$folder.output
 treat<-datatable$treatments
 ntreattemp<-length(treat)
 ntf<-length(expCond$tfNames)
 factorIds<-expCond$factorNames
 nfpexp<-inputparameters$nfpexp
 nFactor<-inputparameters$nFactor
 tfids<-expCond$tfNames
 timevector<-inputparameters$timeVector
 ntime<-inputparameters$ntime
 nplate<-inputparameters$nplate
 TFdis<-inputparameters$TFdis
 tabletemp<-matrix(NA,nrow=(nrepeat*ntreattemp),ncol=ntf)
 tabletempind<-data.frame(matrix(NA,nrow=(nrepeat*ntreattemp),
   ncol=(ntf+nFactor+5))) 
 for (dd in 1:length(treat)){
 	for (cc in 1:ntf) {
		m1 = (datatable$expList[,2] == expCond$tfNames[cc])
		m2 = (datatable$expList[,ncol(datatable$expList)] == treat[dd])
		m3 = m1*m2
		pos = which(as.logical(m1*m2))
		if(length(pos) == 1){
			id_no = datatable$expList[pos,1]
			values= datatable$data[which(datatable$map == as.numeric(id_no))]
			if(length(values)==nrepeat){
				tvalues<-values
			} else{
				tvalues<-c(values,rep(NA,(nrepeat-length(values))))
			}
			tabletemp[(1+(dd-1)*nrepeat):(dd*nrepeat),cc] = as.numeric(tvalues)
			} 
            if (length(pos)>1){
			platei<-ifelse(index%%nplate==0,nplate,index%%nplate)
			id_no = datatable$expList[pos[platei],1]
			values= datatable$data[which(datatable$map == as.numeric(id_no))]
			if(length(values)==nrepeat){
				tvalues<-values
			} else{
				tvalues<-c(values,rep(NA,(nrepeat-length(values))))
			}
			tabletemp[(1+(dd-1)*nrepeat):(dd*nrepeat),cc] = as.numeric(tvalues)
		} 
		if ((length(pos) != 1) && (TFdis!="DIF")){
			print("Warning: Identifier Number Is Not Unique Or Does Not Exist");
			print(c(pos,treat[dd],expCond$tfNames[cc]))
		}
	}
 }
 header<-c(tfids,"repeat","time","plate","treatment","experiment",factorIds)
 colnames(tabletempind)<-header
 treatord<-rep(treat,each=nrepeat)
 factorord<-c()
 if (nFactor>1){
 	for (i in 1:nFactor){
 		tmpfactor<-rep(datatable$levels[,i],each=nrepeat)
      	factorord<-cbind(factorord, as.character(tmpfactor))
 	}
 } else {
		tmpfactor<-rep(datatable$levels,each=nrepeat)
      	factorord<-cbind(factorord, as.character(tmpfactor))
 }
 tabletemp<-apply(tabletemp,c(1,2),function(x){
	r<-ifelse(is.na(x),NA,as.numeric(x))
 })
 tabletempind[,1:ntf]<-tabletemp
 tabletempind[,which(header=="repeat")] =rep(1:nrepeat,ntreattemp)
 tabletempind[,which(header=="time")] = timevector[ifelse(counter%%(ntime)==0,
	(ntime-(counter%%ntime)),counter%%ntime)]
 tabletempind[,which(header=="plate")] = 1
 tabletempind[,which(header=="treatment")] =treatord
 tabletempind[,which(header=="experiment")] = ceiling(counter/nfpexp)
 tabletempind[,(which(header=="experiment")+1):ncol(tabletempind)]<-factorord
 return(tabletempind)
}
 	
# takelog
#=============================================================================================

takelog<-function(datatable=dataAll,expCond=expCond,inputparameters=inputparameters){
  folder.output=inputparameters$folder.output
  ntf<-inputparameters$nTF
  tabletempind<-datatable$tabletemp
  tabletemp<-tabletempind[,1:ntf]
  tabletemp<-apply(tabletemp,c(1,2),function(x){
	r<-ifelse(is.na(x),NA,as.numeric(x))
  })
  ncoltabletemp<-ncol(tabletemp)
  nrowtabletemp<-nrow(tabletemp)
  tabletemp1<-tabletemp
  if (min(tabletemp1,na.rm=TRUE)<=0){
    #Avoid negative number after log	
	tabletemp2<-apply(tabletemp1,c(1,2),function(x){
		r<-ifelse(x==0,0,(x-min(tabletemp1,na.rm=TRUE)+2))
		return(r)
	})
  } else {
	tabletemp2<-tabletemp1
  }
  tablelog<-apply(tabletemp2,c(1,2),function(x){
	r<-ifelse(x==0,0,log(x,2))
	return(r)
  })
  return(tablelog)
}

# backgroundCorr
#=============================================================================================

backgroundCorr<-function(datatable=dataAll,expCond=expCond,inputparameters=inputparameters)
  {
  folder.output=inputparameters$folder.output
  ntf<-inputparameters$nTF
  nplate<-inputparameters$nplate
  ntime<-inputparameters$ntime
  treattemp<-datatable$treatments
  #listtemp<-datatable$expList
  tf_ids<-expCond$tfNames
  tablelog<-datatable$tablelog
  tabletemp<-datatable$tabletemp
  conf.bg<-inputparameters$confBg
  ntreat<-inputparameters$ntreat
  timevector<-inputparameters$timeVector
  #treatment_ids<-expCond$treatmentNames
  summary<-datatable[["Summary"]]
  summaryBcg<-summary
  if (nplate>1){
	 summaryBcg<-summary[which(summary[,2]==1),]
  }
  bg<-inputparameters$bg
  conf.bg<-inputparameters$confBg
  tf.ct<-inputparameters$tfControl
  tf.nodna<-inputparameters$tfBlank
  nFactor<-inputparameters$nFactor
  control_row<-which(tf_ids==tf.ct)
  nodna_row<-which(tf_ids==tf.nodna)
  table1temp<-tablelog
  df <- dim(tablelog)
  h1 <- df[1]
  table_high <- matrix(0,h1, 1)
  for (gg in 1:length(treattemp)){
		nrepeat_NO_DNA<-length(which(!is.na(tablelog[which(tabletemp[,"treatment"]==treattemp[gg]),nodna_row])))
		tmqnorm<-qnorm(conf.bg, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
		conf.bf.rep<-tmqnorm/(nrepeat_NO_DNA^0.5)
		table_high[which(tabletemp[,"treatment"]==treattemp[gg]),1] <-
			mean(tablelog[which(tabletemp[,"treatment"]==treattemp[gg]),nodna_row], na.rm = TRUE) + 
			conf.bf.rep*sd(tablelog[which(tabletemp[,"treatment"]==treattemp[gg]),nodna_row], na.rm = TRUE)
	}
  table_high<-table_high[,rep(1,ntf)]
  table_high<-table1temp-table_high
  table_high<-apply(table_high, c(1,2),function(x){ifelse(!is.na(x) && x<0,TRUE,FALSE)})
  table1temp[table_high]<-0
  for (gg in 1:length(treattemp)){
	for (hh in 1:ntf){
      	rowtreatgg<-which(tabletemp[,"treatment"]==treattemp[gg])
		if (all(is.na(table1temp[rowtreatgg,hh]))){
            	summaryBcg[gg,hh+4]<-NA
  		} else {
			rowtreatggnona<-rowtreatgg[!is.na(tabletemp[rowtreatgg,hh])] 
  			if (all(table1temp[rowtreatggnona,hh]==0)){
            		summaryBcg[gg,hh+4]<-0
  			} else{
				irep<-length(rowtreatggnona)
				aboveBcg<-table1temp[rowtreatggnona,hh]!=0
				summaryBcg[gg,hh+4]<-ifelse (sum(aboveBcg)>(0.5*irep),1,0)			
			}
		}
	}
	rowtreatggnona<-rowtreatgg[!is.na(tabletemp[rowtreatgg,control_row])]
	if (all(table1temp[rowtreatggnona,control_row]==0) || summaryBcg[gg,(control_row+4)]==0){
		table1temp[rowtreatgg,]<-0
            if (nplate==1){
            	summary[gg,ntf+5]<-0
		} else {
			summary[which(summary[,3]==treattemp[gg]),ntf+5]<-0
		}
  	} else{
            if (nplate==1){
            	summary[gg,ntf+5]<-1
		} else {
			summary[which(summary[,3]==treattemp[gg]),ntf+5]<-1
		}
      }
  }
  datatable[["Summary"]]<-summary
  datatable[["SummaryBcg"]]<-summaryBcg
  if (bg==1){
  	datatable[["dataBcgCorr"]]<-table1temp
  } else {
  	datatable[["dataBcgCorr"]]<-tablelog
  }
  table2<-tabletemp
  table2[,1:ntf]<-table_high
  bckTableA<-c()
  bckTableB<-c()
  for (nn in 1:ntf){
	for (oo in 1: nrow(table2)){
      	if (table2[oo,nn]){
			bckTableA<-c(bckTableA,TF=tf_ids[nn])
			bckTableB<-rbind(bckTableB,
				table2[oo,(ntf+1):ncol(table2)]
			)
		}
	}
  }
  bckTable<-data.frame(TF=data.frame(bckTableA),bckTableB)
  if (length(bckTable)!=0){
  	colnames(bckTable)<-c("TF",colnames(table2)[(ntf+1):ncol(table2)])	 
  	bckTable<-bckTable[which(bckTable[,1]!=tf.nodna),]
  }
  datatable[["backgroundTable"]]<-bckTable
  return(datatable)
}

# abovePc
#=============================================================================================

abovePc<-function(datatable=dataAll,expCond=expCond,inputparameters=inputparameters)
  {
    folder.output=inputparameters$folder.output
  tf.pc<-inputparameters$tfPc
  tf.ct<-inputparameters$tfControl
  tf.nodna<-inputparameters$tfBlank
  listtemp<-datatable$expList
  treatmentids<-expCond$treatmentNames
  treattemp<-datatable$treatment
  ntreattemp<-length(treattemp)
  dataBcgCorr<-datatable$dataBcgCorr
  tabletempind<-datatable$tabletemp
  datalog<-datatable$tablelog
  tfids<-expCond$tfNames
  positivecontrol_row<-which(tfids==tf.pc)
  control_row<-which(tfids==tf.pc)
  nodna_row<-which(tfids==tf.nodna)
  conf.pc<-inputparameters$confPc
  aboveBcPc<-inputparameters$aboveBcPc
  summary<-datatable$Summary
  dataPcCorr<-dataBcgCorr
  ntf<-inputparameters$nTF
  if (!is.null(tf.pc)&& (tf.pc %in% unique(listtemp[,2]))){
	tpc_table<-matrix(NA,length(treattemp),1)
	for (ff in 1:ntreattemp){
		PositiveControl_vector <- dataBcgCorr[which(tabletempind[,"treatment"]==treattemp[ff]),positivecontrol_row]
		Pcnona<-PositiveControl_vector[!is.na(PositiveControl_vector )]
		ireppc<-length(Pcnona)
		Control_vector <- log(aboveBcPc,2)+datalog[which(tabletempind[,"treatment"]==treattemp[ff]),nodna_row]
		Ctnona<-Control_vector[!is.na(Control_vector )]
		irepct<-length(Ctnona)
		if (sum(Pcnona!=0)>0.5*ireppc && sum(Ctnona!=0)>0.5*irepct){# to account for NA in the data
			output <- t.test(PositiveControl_vector , Control_vector,alternative = "greater",na.action="na.omit", conf.level = conf.pc) 
			pval <- as.matrix(data.frame(output[3]))
			tpc_table[ff,1] <- pval
		} else{
			tpc_table[ff,1] <- 1
            }
      }
	bin<- (tpc_table<=(1-conf.pc))
	for (ii in 1:ntreattemp){
		if(!is.na(match(treatmentids,treattemp)[ii]) & 
			bin[which(treattemp==treatmentids[ii]),positivecontrol_row]=="TRUE"){
            	if (nplate==1){
            		summary[ii,ntf+6]<-1
			} else {
				summary[which(summary[,3]==treattemp[ii]),ntf+6]<-1
			}
		} else {
			if (nplate==1){
            		summary[ii,ntf+6]<-0
			} else {
				summary[which(summary[,3]==treattemp[ii]),ntf+6]<-0
			}
            }
	}
      if (all(summary[,ntf+6]==0)){
      	dataPcCorr<-apply(dataPcCorr,c(1,2),function(x){x<-NA})
	}
  } else {
	summary[,ntf+6]<-NA
  }
  datatable[["dataPcCorr"]]<-dataPcCorr
  datatable[["Summary"]]<-summary
  return ( datatable)
}

# rmOutlier
#=============================================================================================

rmOutlier<-function(datatable=dataAll,expCond=expCond,inputparameters=inputparameters)
  {
  folder.output=inputparameters$folder.output
  table1<-datatable$dataPcCorr
  treattemp<-datatable$treatments
  listtemp<-datatable$expList
  tf_ids<-expCond$tfNames
  tabletempind<-datatable$tabletemp
  conf.out<-inputparameters$confOutlier
  tf.nodna<-inputparameters$tfBlank
  ntf<-inputparameters$nTF
  if (!all(is.na(table1))){
  	nadetection<-apply(table1,c(1,2),function(x){
		r<-ifelse(is.na(x),0,1)
		return(r)
	})
	if (sum(nadetection)!=0){
		df <- dim(table1)
		h1 <- df[1]
		w1 <- df[2]
		table_low <- matrix(0,h1, w1)
		table_high <- matrix(0,h1, w1)
		table2 <- table1
		for (kk in 1:length(treattemp)){
			for (hh in 1:nTF){
				outdata<-table1[which(tabletempind[,"treatment"]==treattemp[kk]),hh]
				if(tf_ids[hh]%in%listtemp[,2]){	
					tmpmedian<-median(table1[which(tabletempind[,"treatment"]==treattemp[kk]),hh],na.rm=TRUE)
					tmpsd<-sd(table1[which(tabletempind[,"treatment"]==treattemp[kk]),hh],na.rm=TRUE)
					tmqnorm<-qnorm(conf.out, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
					nrepeat_hh<-length(which(!is.na(table1[which(tabletempind[,"treatment"]==treattemp[kk]),hh])))
		            conf.out.rep<-tmqnorm/(nrepeat_hh^0.5)
					table_low[which(tabletempind[,"treatment"]==treattemp[kk]),hh] <- 
						 tmpmedian - conf.out.rep*tmpsd
					table_high[which(tabletempind[,"treatment"]==treattemp[kk]),hh] <- 
						tmpmedian + conf.out.rep*tmpsd
				} else {
					table_low[which(tabletempind[,"treatment"]==treattemp[kk]),hh] <- outdata
					table_high[which(tabletempind[,"treatment"]==treattemp[kk]),hh] <- outdata
				}
			}
        }	
        table_low<-table_low-table2
        table_high<-table2-table_high
        table_low<-apply(table_low,c(1,2),function(x){
				if (is.na(x)){
					r<-NA
				} else {
					r<-ifelse(x>0,TRUE,FALSE)
				}
				return(r)
			})
            table_high<-apply(table_high,c(1,2),function(x){
				if (is.na(x)){
					r<-NA
				} else {
					r<-ifelse(x>0,TRUE,FALSE)
				}
				return(r)
			})
            outlierlist<-apply(table_high+table_low,c(1,2),function(x){
				if (all(is.na(x))){
					r<-NA
				} else {
					r<-sum(x,na.rm=TRUE)
				}
				return(r)
			})
            outlierlist<-apply(outlierlist,c(1,2),function(x){ifelse(x==0 || is.na(x),FALSE,TRUE)})
            table2[outlierlist]<-0
            datatable[["dataRmOutlier"]]<-table2
			outlierTableA<-c()
			outlierTableB<-c()
		if (any(outlierlist)){
  			for (nn in 1:ncol(outlierlist)){
				for (oo in 1: nrow(outlierlist)){
      				if (outlierlist[oo,nn]){
						outlierTableA<-c(outlierTableA,tf_ids[nn])
						outlierTableB<-rbind(outlierTableB,
							tabletempind[oo,(ntf+1):ncol(tabletempind)]
						)
					}
				}
  			}
			outlierTable<-data.frame(TF=data.frame(outlierTableA),outlierTableB)
  			colnames(outlierTable)<-c("TF",colnames(tabletempind)[(ntf+1):ncol(tabletempind)])	 
  			outlierTable<-outlierTable[-which(outlierTable[,1]==tf.nodna),]
			datatable[["outlierTable"]]<-outlierTable 
		} else {
			datatable[["outlierTable"]]<-c()
		}
	} 
  }else {
	table2<-apply(table1,c(1,2),function(x){x<-NA})
	datatable[["outlierTable"]]<-c()
  }
  table2ind<-data.frame(table2,tabletempind[,(nTF+1):ncol(tabletempind)])
  datatable[["dataRmOutlier"]]<-table2
  datatable[["dataRmOutlierind"]]<-table2ind
  return(datatable)
}

# tableGen
#=============================================================================================

tableGen<-function(inputparameters=inputparameters,expCond=expCond,comp=comp){
  folder.output=inputparameters$folder.output
  counter=0
  files=0
  counternorm=1
  map<-inputparameters$condMap
  TFdis<-inputparameters$TFdis
  tfloc<-inputparameters$tfloc
  namefactor<-inputparameters$namesFactors
  expList<-inputparameters$condList
  timevector<-inputparameters$timeVector
  nrepeat<-inputparameters$nrepeat
  nTF<-inputparameters$nTF
  ntreat<-inputparameters$ntreat
  nFactor<-inputparameters$nFactor
  nplate<-inputparameters$nplate
  ntime<-inputparameters$ntime
  nexp<-inputparameters$nexp
  tfids<-expCond$tfNames
  exp2remove<-inputparameters$exp2remove
  if (comp=="Fluc"){
	data<-inputparameters$filesFluc
  } else {
	data<-inputparameters$filesGluc
  }
  for (bb in 1:length(data)){
    arrangeFile<-FALSE
	indSeq<-seq(1,length(data),by=nplate)
	if (nplate>1 && bb%in%indSeq){
	    dataIn<-data[bb:(bb-1+nplate)]
		ordMapIn<-ceiling(bb/(ntime*nplate*nexp))
        mapIn<-map[(1+(ordMapIn-1)*nplate):(ordMapIn*nplate)]
      	expListIn<-expList[(1+(ordMapIn-1)*nplate):(ordMapIn*nplate)]
        files=files+1
		arrangeFile<-TRUE
	} 
	if (nplate==1){
		dataIn<-data[bb]
		ordMapIn<-ceiling(bb/(ntime*nplate))
		ordPlate<-ifelse((bb%%nplate)==0,nplate,bb%%nplate)
        mapIn<-map[(ordPlate+(ordMapIn-1)*nplate)]
		expListIn<-expList[(ordPlate+(ordMapIn-1)*nplate)]
		files=files+1
		arrangeFile<-TRUE
	}
	if (arrangeFile){
      	dataAll<-readData(dataIn=dataIn,mapIn=mapIn,expListIn=expListIn,plateIn=nplate,
  	  		inputparameters)
      	if (tfloc=="FIX"|| (nplate>1 && bb%in%indSeq)){
			if (nFactor==1){
				expListtemp<-cbind(dataAll$expList,dataAll$expList[,3])
                dataAll$expList<-expListtemp
			} else {
				expList2<-dataAll$expList
				expListtemp2<-apply(expList2, 1, function(x){
			  		temp<-as.character(x[3])
			  		for (i in 4:length(x)){
						temp<-paste(temp,as.character(x[i]),sep=".")
			  		}
			  		return(temp)
            		})
	      		expListtemp<-cbind(expList2,expListtemp2)
				dataAll$expList<-expListtemp
      		}
            maptempIds<-unique(as.vector(dataAll$map))
			treattemp<-unique(expListtemp[na.omit(match(maptempIds,as.numeric(expListtemp[,1]))),
				ncol(expListtemp)])
			leveltemp<-expListtemp[match(treattemp,expListtemp[,ncol(expListtemp)]),3:(2+length(namefactor))]
      		dataAll[["treatments"]]<-treattemp	
            	dataAll[["levels"]]<-leveltemp	
			summary<-fillSummary(counter=files,inputparameters=inputparameters,datatable=dataAll,
 				expCond=expCond,expListtemp=expListIn)
			dataAll[["Summary"]]<-summary 
      		dataOrder<-orderTables(datatable=dataAll,expCond=expCond,inputparameters=inputparameters,
				counter=files,index=bb)	
      		dataAll[["tabletemp"]]<-dataOrder
			tableLog<-takelog(datatable=dataAll,expCond,inputparameters)
			dataAll[["tablelog"]]<-tableLog
			dataAll<-backgroundCorr(datatable=dataAll,expCond=expCond,inputparameters=inputparameters)
            	dataAll<-abovePc(datatable=dataAll,expCond=expCond,inputparameters=inputparameters)	
			dataAll<-rmOutlier(datatable=dataAll,expCond=expCond,inputparameters=inputparameters)
		}
      	if (bb==1){
            finalTable<-dataAll$tabletemp
            bcgCorrected<-dataAll$dataBcgCorr
			prenormTable<-dataAll$dataRmOutlierind
      		summaryTable<-dataAll$Summary
      		summaryBcgTable<-dataAll$SummaryBcg
            backgroundList<-dataAll$backgroundTable
            outlierList<-dataAll$outlierTable
      	} else {
 			finalTable<-rbind(finalTable,dataAll$tabletemp)
 			bcgCorrected<-rbind(bcgCorrected,dataAll$dataBcgCorr)
     		prenormTable<-rbind(prenormTable,dataAll$dataRmOutlierind)
			summaryTable<-rbind(summaryTable,dataAll$Summary)
			summaryBcgTable<-rbind(summaryBcgTable,dataAll$SummaryBcg)
			backgroundList<-rbind(backgroundList,dataAll$backgroundTable)
            outlierList<-rbind(outlierList,dataAll$outlierTable)
		}
  	}
  }
  header<-c(tfids,"repeat","time","plate","treatment","experiment",namefactor)
  headersummaryTable<-c("Experiment","Plate","Treatment","Time",
  paste(tfids,"in_experiment",sep="_"),"TA_above_bcg","PC_above_bcg")
  headersummaryBcgTable<-c("Experiment","Plate","Treatment","Time",
  paste(tfids,"aboveBcg",sep="_"))
  colnames(finalTable)<-header
  colnames(prenormTable)<-header
  colnames(summaryTable)<-headersummaryTable
  summaryBcgTable<-summaryBcgTable[,1:(ncol(summaryBcgTable)-2)]
  colnames(summaryBcgTable)<-headersummaryBcgTable
  finalTable<-finalTable[
	order(as.numeric(finalTable[,"experiment"]),
	as.numeric(finalTable[,"time"]),
	finalTable[,"treatment"],
	as.numeric(finalTable[,"repeat"])),
  ]
  prenormTable<-prenormTable[
	order(as.numeric(prenormTable[,"experiment"]),
	as.numeric(prenormTable[,"time"]),
	prenormTable[,"treatment"],
	as.numeric(prenormTable[,"repeat"])),
  ]
  prenormTabletmp <- prenormTable
  prenormTabletmp[,1:nTF]<-apply(prenormTable[,1:nTF],c(1,2),function(x){
	if (is.na(x) || x==0){
		r=0
	} else {
		r<-2^x
	} 
	return(r)
  })
  prenormTabletmp[is.na(prenormTable)]<-NA
  prenormTable<-prenormTabletmp  
  AllTables<-list(finalTable=finalTable,prenormTable=prenormTable,summaryTable=summaryTable,
	summaryBcgTable=summaryBcgTable,backgroundList=backgroundList,outlierList=outlierList)
  write.table(finalTable, file = paste(folder.output,"/RawData.txt",sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = NA, qmethod = c("escape", "double"))
  write.table(backgroundList, file = paste(folder.output,"/DataBelowBackground.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = NA, qmethod = c("escape", "double"))
  write.table(outlierList, file = paste(folder.output,"/Outliers.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = NA, qmethod = c("escape", "double"))
  write.table(prenormTable, file =paste(folder.output,"/TransformData.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = NA, qmethod = c("escape", "double"))
  return(AllTables)
}

# withInArrayShifting
#=========================================================================================

withInArrayShifting<-function(datatable=dataAll,expCond=expCond,inputparameters=inputparameters)
  {
  folder.output=inputparameters$folder.output
  prenormind<-datatable$prenormTable
  nTF<-inputparameters$nTF
  treat<-expCond$treatmentNames
  prenorm<-prenormind[,1:nTF]
  tmpnorm<- prenorm
  ntreat<-length(treat)
  timevector<-inputparameters$timeVector
  nexp<-inputparameters$nexp
  tWithin<-inputparameters$tWithin
  if (tWithin && !all(is.na(prenorm))){	
		for (ll in 1:nexp){
  			for (mm in 1:ntreat){
				alldatai<-prenormind[which(
					prenormind[,"experiment"]==ll &
					prenormind[,"treatment"]==treat[mm]),
				]
				for (kk in 1:(nTF)){
					timesum<-tapply(alldatai[,kk],alldatai[,"time"],function(x){
							sum(x,na.rm=TRUE)
					})
                    timeshift<-timevector[which(timesum!=0)][1]
					#JIA: get the correct median of the section.
					TF_repeats <- alldatai[which(alldatai!=0 & alldatai[,"time"]==timeshift & col(alldatai) == kk, arr.ind=TRUE)]
					medianalli<-median(as.numeric(TF_repeats),na.rm=TRUE)
                    if (is.na(timeshift)){
						tmpnorm[which(prenormind[,"experiment"]==ll &
						prenormind[,"treatment"]==treat[mm]),kk]<-0
					} else {
						tmpnorm[which(prenormind[,"experiment"]==ll &
							prenormind[,"treatment"]==treat[mm]),kk]<-
							alldatai[,kk]*medianalli/alldatai[which(alldatai[,"time"]==timeshift),kk
						]
					}
				}
			}
		}
  }
  tmpnorm<-apply(tmpnorm,c(1,2),function(x){
	r<-ifelse(x=="NaN" || x=="Inf",0,x)
	return(r)
  })
  tmpnorm[is.na(prenormind[,1:nTF])]<-NA
  dataWithInShift<-data.frame(tmpnorm,prenormind[,(nTF+1):ncol(prenormind)])
  colnames(dataWithInShift)<-colnames(prenormind)
  datatable[["dataWithInShift"]]<-dataWithInShift
  return(datatable)
}

# Glucnorm
#============================================================================================

Glucnorm<-function(dataGluc=dataAllGluc,dataFluc=dataAllFluc,
	inputparameters=inputParameters,expCond=expCond){
	  folder.output=inputparameters$folder.output
  if (inputparameters$tGluc){
	  datagluc<-dataGluc$dataWithInShift
      datafluc<-dataFluc$dataWithInShift
      dataflucgluc<-datafluc
      nTF<-inputparameters$nTF
      ntreat<-inputparameters$ntreat
      treat<-expCond$treatmentNames
      timevector<-inputparameters$timeVector
      nexp<-inputparameters$nexp
      ntime<-inputparameters$ntime
	for (aa in 1:nexp){
		for (bb in 1:ntreat){
			Gluc0<-datagluc[which(
				datagluc[,"experiment"]==aa &
				datagluc[,"treatment"]==treat[bb] &
				datagluc[,"time"]==timevector[1]),1:nTF
			]
			Gluc0<-Gluc0[rep(1:nrow(Gluc0),ntime),]
	      	dataflucgluc[which(
				dataflucgluc[,"experiment"]==aa &
				dataflucgluc[,"treatment"]==treat[bb]),1:nTF
			]<-datafluc[which(
				datafluc[,"experiment"]==aa &
				datafluc[,"treatment"]==treat[bb]),1:nTF
			]/Gluc0
		}
	}
  	dataflucgluc<-apply(dataflucgluc[,1:nTF],c(1,2),function(x){
		r<-ifelse(x=="NaN" || x==Inf,0,x)
		return(r)
	})
	dataflucgluc[is.na(datafluc[,1:nTF])]<-NA
	dataflucgluc<-data.frame(dataflucgluc,datafluc[,(nTF+1):ncol(datafluc)])
  	colnames(dataflucgluc)<-colnames(datafluc)
  } else {
	dataflucgluc<-dataFluc$dataWithInShift
  }	
  dataFluc[["dataFlucGluc"]]<-dataflucgluc
  return(dataFluc)
}

# normTA
#=============================================================================================

normTA<-function(datatable=dataFluc,expCond=expCond,inputparameters=inputparameters)
  {
    folder.output=inputparameters$folder.output
  dataWithInInd<-datatable$dataFlucGluc
  nTF<-inputparameters$nTF
  treat<-expCond$treatmentNames
  dataWithIn<-dataWithInInd[,1:nTF]
  datanormTA<- dataWithIn
  ntreat<-length(treat)
  tf.ct<-inputparameters$tfControl
  tfids<-expCond$tfNames
  control_row<-which(tfids==tf.ct)
  timevector<-inputparameters$timeVector
  ntimes<-length(timevector)
  nexp<-inputparameters$nexp
  nrepeats<-inputparameters$nrepeat
  if (!all(is.na(dataWithIn))){
	for (ff in 1:nTF){	
		for (gg in 1:nexp){
  			for (hh in 1:ntreat){
				for (ii in 1:ntimes){
					dataTAi<-dataWithInInd[which(
							dataWithInInd[,"experiment"]==gg &
							dataWithInInd[,"treatment"]==treat[hh] &
							dataWithInInd[,"time"]==timevector[ii]),
							control_row
					]
                    if (all(dataTAi[!is.na(dataTAi)]==0)|| all (is.na(dataTAi))){
						medianTAi<-0
					} else {
						medianTAi<-median(dataTAi[which(dataTAi!=0)],
						na.rm=TRUE)
					}
					tmp<-tryCatch(
						dataWithInInd[which(
							dataWithInInd[,"experiment"]==gg &
							dataWithInInd[,"treatment"]==treat[hh] &
							dataWithInInd[,"time"]==timevector[ii]),
							ff
						]/medianTAi,
						error=function(e) e
					)
					if(!inherits(tmp, "error")){
    						datanormTA[which(
							dataWithInInd[,"experiment"]==gg &
							dataWithInInd[,"treatment"]==treat[hh] &
							dataWithInInd[,"time"]==timevector[ii]),
							ff]<-tmp
 					} else {
						 datanormTA[which(
							dataWithInInd[,"experiment"]==gg &
							dataWithInInd[,"treatment"]==treat[hh] &
							dataWithInInd[,"time"]==timevector[ii]),
							ff]<-0
					}	
				}
			}
		}
	}
  }
  datanormTA[is.na(dataWithInInd[,1:nTF])]<-NA
  datanormTAind<-data.frame(datanormTA,dataWithInInd[,(nTF+1):ncol(dataWithInInd)])
  colnames(datanormTAind)<-colnames(dataWithInInd)
  datatable[["dataTAnormind"]]<-datanormTAind
  write.table(datanormTAind, file = paste(folder.output,"/NormalizedTA.txt",sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = NA, qmethod = c("escape", "double"))
  return(datatable)
}

# withInAndTA
#============================================================================================	
			
withInAndTA<-function(dataGluc=dataAllGluc,dataFluc=dataAllFluc,
	inputparameters=inputparameters,expCond=expCond){
  folder.output=inputparameters$folder.output
  if (!is.na(filesGluc)[1]){
  	dataGluc<-withInArrayShifting(datatable=dataGluc,
		expCond=expCond,inputparameters=inputparameters)
  	plotWithInArrayShiftingData(dataAll=dataGluc,inputparameters=inputParameters,
		expCond=expCond,comp="Gluc")
  }
  dataFluc<-withInArrayShifting(datatable=dataFluc,
	expCond=expCond,inputparameters=inputparameters)
  plotWithInArrayShiftingData(dataAll=dataFluc,inputparameters=inputParameters,
	expCond=expCond,comp="Fluc")
  dataFluc<-Glucnorm(dataGluc=dataGluc,dataFluc=dataFluc,
	inputparameters=inputParameters,expCond=expCond)
  dataFluc<-normTA(datatable=dataFluc,
	expCond=expCond,inputparameters=inputparameters)
  plotnormTA(dataAll=dataFluc,inputparameters=inputParameters,
	expCond=expCond) 
  library(plyr)
 # In case of control treatment presents
  nameUntreated<-inputparameters$nameUntreated
  if (any(nameUntreated!="none")){
    # required inputs
	normTA<-dataFluc$dataTAnormind
	nrepeats<-inputparameters$nrepeat
	ntreatments<-inputparameters$ntreat
	nFactor<-inputparameters$nFactor
	treatmentNames<-expCond$treatmentNames
	if (nFactor==1){
		normTAvsCtr<-ddply(normTA[which(normTA[,"treatment"]==nameUntreated),],.(experiment,time),function(x){
			apply(x[,1:nTF],2,function(y){
				average<-mean(y,na.rm=TRUE)
				r<-ifelse(average=="NaN",NA,average)
				return(r)
				})
		})
		normTAvsCtr<-normTAvsCtr[rep(seq(1,nrow(normTAvsCtr),by=1),each=(nrepeats*ntreatments)),]
		normTAvsCtrl<-normTA[,1:nTF]/normTAvsCtr[,3:ncol(normTAvsCtr)]
		normTAvsCtrl<-apply(normTAvsCtrl,c(1,2),function(x){
			r<-ifelse((x=="NaN" || x=="inf"),NA,x)
			return(r)
		})
	} else {
		# identify the factors that have a control
		fid<-which(nameUntreated!="none")
		# generate a temporal table with the corresponding control ids
		treattemp<-NULL
		for (rr in 1:nrow(normTA)){
			treattemp<-rbind(treattemp,
				unlist(strsplit(as.character(normTA$treatment[rr]),"\\.")))
		}
		factortemp<-data.frame(treattemp)
		for (tt in 1:length(fid)){
			treattemp[,fid[tt]]<-nameUntreated[fid]
		}
		treattmp<-apply(treattemp,1,function(x){paste(x,collapse=".")})
		normTA[,"treattmp"]<-treattmp
		controlNames<-unique(normTA[,"treattmp"])
		# identify the location of the temporal treatments
		postreattmp<-grep("treattmp",colnames(normTA))
		locUnt<-apply(normTA,1,function(x){
			id<-x[postreattmp]
			rr<-ifelse(id%in%controlNames,TRUE,FALSE)
			return(rr)
		})
		# calculate the averages of the control
		normTAvsCtr<-ddply(normTA[locUnt,],.(experiment,treatment,time),function(x){
			apply(x[,1:nTF],2,function(y){
				average<-mean(y,na.rm=TRUE)
				r<-ifelse(average=="NaN",NA,average)
				return(r)
				})
		})
		# Normalize by the control data
		normTAvsCtrl<-ddply(normTA,.(experiment,treatment,time),function(x){
			ctrtmp<-normTAvsCtr[which(normTAvsCtr$experiment==x$experiment[1] & 
				  normTAvsCtr$time==x$time[1] &
				  normTAvsCtr$treatment==x$treattmp[1]),]
			zz<-x[,1:nTF]/ctrtmp[,match(colnames(normTA,normTAvsCtr))]
			res<-apply(zz,c(1,2),function(x){
				r<-ifelse((x=="NaN" || x=="inf"),NA,x)
				return(r)
			})
			return(res)
		})
	}
	normTAvsCtrl<-data.frame(normTAvsCtrl,normTA[,(nTF+1):ncol(normTA)])
	dataFluc[["normTAvsCtrl"]]<-normTAvsCtrl
	plotnormTAvsControl(dataAll=dataFluc,inputparameters=inputParameters,
		expCond=expCond)
  }
  return(dataFluc)
}

# shiftall
#============================================================================================	

 shiftall<-function(dataInd=dataShiftInd,data=dataShift,ntimes=ntimes,
			nexp=nexp,timevector=timevector,nTF=nTF){
	datai<-data[which(
		dataInd[,"time"]==timevector[1]),
	]
	medianShifti<-apply(datai,2,function(x){
		median(x[which(x!=0)],na.rm=TRUE)
	})
	for (ff in 1:nexp){
		den<-dataInd[which(
				dataInd[,"experiment"]==ff &
				dataInd[,"time"]==timevector[1]),
		1:nTF]
		den<-den[rep(1:nrow(den),ntimes),]
		tmp<-dataInd[which(
			dataInd[,"experiment"]==ff),
			1:nTF]/den
		data[which(
			dataInd[,"experiment"]==ff),]<-
			t(apply(tmp,1,function(x){x*medianShifti})
		)    		
	}
	if (any(medianShifti[!is.na(medianShifti)]==0)){
	    tfid00<-which(medianShifti==0)
		for (mm in 1:length(tfid00)){
			timesum<-tapply(data[,tfid00],data[,"time"],function(x){
				sum(x,na.rm=TRUE)})
            timeshift<-timevector[which(timesum!=0)][1]
			seldata<-alldatai[which(alldatai[,"time"]==timeshift),kk]
			medianShiftmm<-median(seldata[which(seldata!=0)],na.rm=TRUE)
			for (ff in 1:nexp){
				den<-dataInd[which(
					dataInd[,"experiment"]==ff &
					dataInd[,"time"]==timeshift),tfid00[mm]]
				den<-den[rep(1:nrow(den),ntimes),]
				tmp<-dataInd[which(
					dataInd[,"experiment"]==ff),tfid00[mm]]/den
				data[which(
					dataInd[,"experiment"]==ff),tfid00[mm]]<-
						tmp*medianShiftmm
			}   
		}
	}
    data<-data.frame(data,dataInd[,(nTF+1):ncol(dataInd)])
	colnames(data)<-colnames(dataInd)
	return(data)
 }

# shiftpartial
#============================================================================================	

 shiftpartial<-function(dataInd=dataShiftInd,data=dataShift,
	nexp=nexp,treat=treat,ntreat=ntreat,ntimes=ntimes,
	timevector=timevector,nTF=nTF){
 for (gg in 1:ntreat){
	datai<-data[which(
		dataInd[,"time"]==timevector[1] &
		dataInd[,"treatment"]==treat[gg]),
	]
	medianShifti<-apply(datai,2,function(x){
		median(x[which(x!=0)],na.rm=TRUE)
	})
	for (ff in 1:nexp){
		den<-dataInd[which(
			dataInd[,"experiment"]==ff &
			dataInd[,"time"]==timevector[1]&
			dataInd[,"treatment"]==treat[gg]
		),1:nTF]
		den<-den[rep(1:nrow(den),ntimes),]
	    	tmp<-dataInd[which(
			dataInd[,"experiment"]==ff &
			dataInd[,"treatment"]==treat[gg]
		),1:nTF]/den
		data[which(
			dataInd[,"experiment"]==ff &
			dataInd[,"treatment"]==treat[gg]
		),]<-t(apply(tmp,1,function(x){x*medianShifti}))
	}
	if (any(medianShifti[!is.na(medianShifti)]==0)){
	    tfid00<-which(medianShifti==0)
		for (mm in 1:length(tfid00)){
			timesum<-tapply(data[,tfid00],data[,"time"],function(x){
				sum(x,na.rm=TRUE)})
            timeshift<-timevector[which(timesum!=0)][1]
			seldata<-alldatai[which(alldatai[,"time"]==timeshift),kk]
			medianShiftmm<-median(seldata[which(seldata!=0)],na.rm=TRUE)
			for (ff in 1:nexp){
				den<-dataInd[which(
					dataInd[,"experiment"]==ff &
					dataInd[,"time"]==timeshift&
					dataInd[,"treatment"]==treat[gg]
				),tfid00[mm]]
				den<-den[rep(1:nrow(den),ntimes),]
				tmp<-dataInd[which(
					dataInd[,"experiment"]==ff &
					dataInd[,"treatment"]==treat[gg]
					),tfid00[mm]]/den
				data[which(
					dataInd[,"experiment"]==ff &
					dataInd[,"treatment"]==treat[gg]
				),tfid00[mm]]<-tmp*medianShiftmm
			}
		}
	}
  }
  data<-apply(data,c(1,2),function(x){
		r<-ifelse(x=="NaN" || x=="Inf",0,x)
		return(r)
  })
  data<-data.frame(data,dataInd[,(nTF+1):ncol(dataInd)])
	colnames(data)<-colnames(dataInd)
  return(data)
 }

# normBetween
#============================================================================================
normBetween<-function(dataFluc=dataAllFluc,
	inputparameters=inputParameters,expCond=expCond){
 folder.output=inputparameters$folder.output
 datat0NormInd<-dataFluc$dataTAnormind
 nTF<-inputparameters$nTF
 treat<-expCond$treatmentNames
 datat0Norm<-datat0NormInd[,1:nTF]
 ntreat<-inputparameters$ntreat
 timevector<-inputparameters$timeVector
 ntimes<-length(timevector)
 nexp<-inputparameters$nexp
 nrepeats<-inputparameters$nrepeat
 nfactor<-inputparameters$nFactor
 namesFactors<-inputParameters$namesFactors
 levelsIds<-expCond$levelsNames
 tIn<-inputparameters$tIn
 expList<-expCond$expList
 if (inputparameters$t0Norm && !all(is.na(datat0Norm))){	
	for (gg in 1:nexp){
  		for (hh in 1:ntreat){
		treatment_data <- datat0NormInd[which(
			datat0NormInd[,"experiment"]==gg & 
			datat0NormInd[,"treatment"]==treat[hh]), 1:nTF]
		t0_treatment_data <- datat0NormInd[which(
			datat0NormInd[,"experiment"]==gg & datat0NormInd[,"treatment"]==treat[hh] & 
			datat0NormInd[,"time"]==timevector[1]), 1:nTF]
		size_factor <- nrow(treatment_data)/nrow(t0_treatment_data)
		sized_t0_treatment_data<-do.call("rbind", replicate(size_factor, 
			t0_treatment_data, simplify = FALSE))
		datat0Norm[which(datat0NormInd[,"experiment"]==gg & datat0NormInd[,"treatment"]==treat[hh]),1:nTF] <-
			treatment_data / sized_t0_treatment_data
		}
	}
	datat0Norm<-apply(datat0Norm,c(1,2),function(x){
		r<-ifelse(x=="NaN" || x=="Inf",0,x)
		return(r)
	})
	datat0Norm[is.na(datat0NormInd[,1:nTF])]<-NA
  } else {
	datat0Norm<-datat0NormInd[,1:nTF]
  }	
  datat0Norm<-data.frame(datat0Norm,datat0NormInd[,(nTF+1):ncol(datat0NormInd)])
  colnames(datat0Norm)<-colnames(datat0NormInd)
  dataFluc[["datat0Norm"]]<-datat0Norm
  dataShiftInd<-datat0Norm
  dataShift<-dataShiftInd[,1:nTF]
  levelsNames<-expCond$levelsNames
  if (tIn!=0 && !all(is.na(dataShift))){
	if (all(tIn==1)){
		dataShiftInd<-shiftall(dataInd=dataShiftInd,data=dataShift,ntimes=ntimes,
			nexp=nexp,timevector=timevector,nTF=nTF)
    } else {
		fid<-which(tIn==2)
		temptreatid<-match(namesFactors[fid],colnames(dataShiftInd))
        if ((nfactor>1) && (length(temptreatid)>1)){
				treattmp<-apply(dataShiftInd[,temptreatid],1,function(x){
					paste0(x,collapse=".")}
				)
        } else {
			treattmp<-dataShiftInd[,temptreatid]
		}
		dataShiftIndtmp<-dataShiftInd
		dataShiftInd0<-dataShiftInd
		dataShiftIndtmp[,"treatment"]<-treattmp
		treat<-unique(treattmp)
		ntreat<-length(treat)	
		dataShiftInd<-shiftpartial(dataInd=dataShiftIndtmp,data=dataShift,
			nexp=nexp,treat=treat,ntreat=ntreat,ntimes=ntimes,
			timevector=timevector,nTF=nTF)
            dataShiftInd[,"treatment"]<-dataShiftInd0[,"treatment"]
	}
	dataShiftInd[,1:nTF]<-apply(dataShiftInd[,1:nTF],c(1,2),function(x){
		r<-ifelse(x=="NaN" || x=="Inf",0,x)
		return(r)
	})
  } else {
	dataShiftInd<-datat0Norm
  }
  dataFluc[["dataShift"]]<-dataShiftInd
  plotnorm(dataAll=dataShiftInd,inputparameters=inputParameters,
	expCond=expCond,pdfname="Shift.pdf")
  #Scale Free
  dataScaleFreeInd<-dataShiftInd
  dataScaleFreeInd[,1:nTF]<-apply(dataScaleFreeInd[,1:nTF],c(1,2),function(x){
	if (x==0 || is.na(x)){
		r<-ifelse(x==0,0,NA)
	} else {
		r<-log(x,2)
	}
  }) 
  dataScaleFree<-dataScaleFreeInd[,1:nTF]
  if ((inputparameters$scaleFree ||inputparameters$scaleMean) && !all(is.na(dataScaleFree))){
	for (cc in 1:nexp){
		datascale<-dataScaleFree[which(
			dataScaleFreeInd[,"experiment"]==cc),
		]
		med.att <- apply(datascale, 2, median,na.rm=TRUE)
		datascale<-sweep(datascale, 2, med.att)
		sdscale<-apply(datascale,2,function(x){
				sd(x,na.rm=TRUE)
		})
		if(!inputparameters$scaleMean){
			datascale<-t(apply(datascale,1,function(x){
				x/sdscale
			}))
		}
		dataScaleFree[which(dataScaleFreeInd[,"experiment"]==cc),]<-datascale
	}
	dataScaleFree<-apply(dataScaleFree,c(1,2),function(x){
		r<-ifelse(x=="NaN" || x=="Inf",0,x)
		return(r)
	})
  } else {
	dataScaleFree<-dataScaleFreeInd[,1:nTF]
  }	
  dataScaleFree<-data.frame(dataScaleFree,dataScaleFreeInd[,(nTF+1):ncol(dataScaleFreeInd)])
  colnames(dataScaleFree)<-colnames(dataScaleFreeInd)
  dataFluc[["dataScaleFree"]]<-dataScaleFree
  write.table(dataScaleFree, file = paste(folder.output,"/AfterAllNormalization.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = NA, qmethod = c("escape", "double"))
  # In case of control treatment present
  nameUntreated<-inputparameters$nameUntreated
    library(plyr)
  if (any(nameUntreated!="none")){
    # required inputs
	normalized<-dataFluc[["dataScaleFree"]]
	nrepeats<-inputparameters$nrepeat
	ntreatments<-inputparameters$ntreat
	nFactor<-inputparameters$nFactor
	treatmentNames<-expCond$treatmentNames
	if (nFactor==1){
		normalizedvsCtr<-ddply(normalized[which(normalized[,"treatment"]==nameUntreated),],.(experiment,time),
			function(x){
				apply(x[,1:nTF],2,function(y){
					tt<-y[y!=0]
					if(length(tt)>0){
						average<-mean(tt,na.rm=TRUE)
						r<-ifelse(average=="NaN",NA,average)
					} else {
						r<-NA
					}
					return(r)
				})
		})
		normalizedvsCtr<-normalizedvsCtr[rep(seq(1,nrow(normalizedvsCtr),by=1),each=(nrepeats*ntreatments)),]
		normalizedvsCtr<-normalized[,1:nTF]-normalizedvsCtr[,3:ncol(normalizedvsCtr)]
		normalizedvsCtr<-apply(normalizedvsCtr,c(1,2),function(x){
			r<-ifelse((x=="NaN" || x=="inf"),NA,x)
			return(r)
		})
	} else {
		# identify the factors that have a control
		fid<-which(nameUntreated!="none")
		# generate a temporal table with the corresponding control ids
		treattemp<-NULL
		for (rr in 1:nrow(normalized)){
			treattemp<-rbind(treattemp,
				unlist(strsplit(normalized$treatment[rr],"\\.")))
		}
		factortemp<-data.frame(treattemp)
		for (tt in 1:length(fid)){
			treattemp[,fid[tt]]<-nameUntreated[fid]
		}
		treattmp<-apply(treattemp,1,function(x){paste(x,collapse=".")})
		normalized[,"treattmp"]<-treattmp
		# identify the location of the temporal treatments
		locUnt<-NULL
		for (ii in 1:length(fid)){
			locUnttmp<-grep(nameUntreated[fid[ii]],normalized$treatment)
			if(ii==1){
				locUnt<-locUnttmp
			} else {
				locUnt<-locUnt[na.omit(match(locUnttmp,locUnt))]
			}
		}
		# calculate the averages of the control
		normalizedvsCtr<-ddply(normalized[locUnt,],.(experiment,treatment,time),function(x){
			apply(x[,1:nTF],2,function(y){
				tt<-y[y!=0]
				if (length(tt)>0){
					average<-mean(y,na.rm=TRUE)
					r<-ifelse(average=="NaN",NA,average)
				} else {
					r<-NA
				}
				return(r)
				})
		})
		# Normalize by the control data
		normalizedvsCtr<-ddply(normalized,.(experiment,treatment,time),function(x){
			ctrtmp<-normalizedvsCtr[which(normalizedvsCtr$experiment==x$experiment[1] & 
				  normalizedvsCtr$time==x$time[1] &
				  normalizedvsCtr$treatment==x$treattmp[1]),]
			zz<-x[,1:nTF]-ctrtmp[,match(colnames(normalized,normalizedvsCtr))]
			res<-apply(zz,c(1,2),function(x){
				r<-ifelse((x=="NaN" || x=="inf"),NA,x)
				return(r)
			})
			return(res)
		})
	}
	normalizedvsCtr<-data.frame(normalizedvsCtr,normalized[,(nTF+1):ncol(normalized)])
	dataFluc[["normalizedvsCtrl"]]<-normalizedvsCtr
	plotnormTAvsControl(dataAll=dataFluc,inputparameters=inputParameters,
		expCond=expCond)
  }
 return(dataFluc)
}

# plotRawDataPDF
#============================================================================================

plotRawDataPDF<-function(dataAll=dataAll,inputparameters=inputParameters,
	expCond=expCond,comp=comp){
  folder.output=inputparameters$folder.output
  valueColors<-inputparameters$valueColors
  library(ggplot2)
  library(reshape)
  library (RColorBrewer)
  rawData<-dataAll$finalTable
  if (comp=="Gluc"){
	pdfname<-"RawDataGluc.pdf"
  } else {
	pdfname<-"RawDataFluc.pdf"
  }
  pdf(file = paste(folder.output,"/",pdfname,sep=""),paper="a4r",width=11,height=8)
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  treatmentids<-expCond$treatmentNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  nFactors<-inputparameters$nFactor
  namesFactors<-inputparameters$namesFactors
  levFactor<-inputparameters$levFactor
  rawData[,"group"]<-rep(rep(seq(1,nrepeat*ntreatment,by=1),ntime),nexp)
  # In linear scale
  if (nFactor==1){
	rawDatatmp<-rawData
	for(i in 1:nTF){
		rawDatatmp[,"TF"]<-rawData[,i]
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		pp<-ggplot(rawDatatmp,aes(x=time,y=TF,
			colour=treatment,group=group))+
			geom_line(size=1)+
			facet_wrap( ~ experiment, ncol=nCol)+
			xlab("Time") +
			ylab("Raw Intensity") +
			labs(title =titlei)+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		if(is.null(valueColors)){
			pp<-pp+scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")
		} else {
			pp<-pp+scale_color_manual(name = "Treatment", values = valueColors)
		}
		print(pp)
	}
  } else {
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(rawData))
		for (kk in 1:levFactor[jj]){
				rawDatatmp<-rawData[which(rawData[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
				for(i in 1:nTF){
					rawDatatmp[,"TF"]<-rawDatatmp[,i]
					nRow=ceiling(-1+sqrt(2+nexp))
					nCol=2+nRow
					titlei<-tf_ids[i]
					pp<-ggplot(rawDatatmp,aes(x=time,y=TF,
						colour=treatment,group=group))+
						geom_line(size=1)+
						facet_wrap( ~ experiment, ncol=nCol)+
						xlab("Time") +
						ylab("Raw Intensity") +
						labs(title =titlei)+
						theme(panel.background = element_rect(fill="transparent"))+
						theme(strip.background=element_rect(fill="white"))+
						theme(axis.ticks=element_line(colour="black"))+
						theme(axis.title=element_text(size=16))+
						theme(panel.border=element_rect(colour="black",fill="transparent"))
					if(is.null(valueColors)){
						pp<-pp+scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")
					} else {
						pp<-pp+scale_color_manual(name = "Treatment", values = valueColors)
					}
					print(pp)
				}
			}	
	}
  }	
  dev.off()
  # In log scale
  pdf(file = paste(folder.output,"/","LogScale",pdfname,sep=""),paper="a4r",width=11,height=8)
  rawDatalog<-data.frame(
	apply(rawData[,1:nTF],c(1,2),function(x){
		r<-ifelse((!is.na(x)||x!=0),log(x,2),x)}),
		rawData[,(nTF+1):ncol(rawData)])
  if (nFactor==1){
		for(i in 1:nTF){
			rawDatatmp[,"TF"]<-rawDatalog[,i]
			nRow=ceiling(-1+sqrt(2+nexp))
			nCol=2+nRow
			titlei<-tf_ids[i]
			pp<-ggplot(rawDatatmp,aes(x=time,y=TF,
				colour=treatment,group=group))+
				geom_line(size=1)+
				facet_wrap( ~ experiment, ncol=nCol)+
				xlab("Time") +
				ylab("Raw Intensity (log scale)") +
				labs(title =titlei)+
				theme(panel.background = element_rect(fill="transparent"))+
				theme(strip.background=element_rect(fill="white"))+
				theme(axis.ticks=element_line(colour="black"))+
				theme(axis.title=element_text(size=16))+
				theme(panel.border=element_rect(colour="black",fill="transparent"))
			if(is.null(valueColors)){
				pp<-pp+scale_colour_brewer(type="qual",palette=1,
					name = "Treatment")
			} else {
				pp<-pp+scale_color_manual(name = "Treatment", values = valueColors)
			}
			print(pp)
		}
  } else {
		for (jj in 1:nFactors){
			posjj<-grep(namesFactors[jj],colnames(rawData))
			for (kk in 1:levFactor[jj]){
				rawDatatmp<-rawDatalog[which(
					rawDatalog[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
				for(i in 1:nTF){
					rawDatatmp[,"TF"]<-rawDatatmp[,i]
					nRow=ceiling(-1+sqrt(2+nexp))
					nCol=2+nRow
					titlei<-tf_ids[i]
					pp<-ggplot(rawDatatmp,aes(x=time,y=TF,
						colour=treatment,group=group))+
						geom_line(size=1)+
						facet_wrap( ~ experiment, ncol=nCol)+
						xlab("Time") +
						ylab("Raw Intensity (log scale)") +
						labs(title =titlei)+
						theme(panel.background = element_rect(fill="transparent"))+
						theme(strip.background=element_rect(fill="white"))+
						theme(axis.ticks=element_line(colour="black"))+
						theme(axis.title=element_text(size=16))+
						theme(panel.border=element_rect(colour="black",fill="transparent"))
					if(is.null(valueColors)){
						pp<-pp+scale_colour_brewer(type="qual",palette=1,
							name = "Treatment")
					} else {
						pp<-pp+scale_color_manual(name = "Treatment", values = valueColors)
					}
					print(pp)
				}
			}	
		}
    }			
  dev.off()
  }

# plotpreNormData
#============================================================================================

plotpreNormData<-function(dataAll=dataAll,inputparameters=inputParameters,
	expCond=expCond,comp=comp){
  folder.output=inputparameters$folder.output
  library(ggplot2)
  library(reshape)
  library (RColorBrewer)
    library(plyr)
  prenormTable<-dataAll$prenormTable
  if (comp=="Gluc"){
	pdfname<-"TransBckOutlierDataGluc.pdf"
  } else {
	pdfname<-"TransBckOutlierDataFluc.pdf"
  }
  pdf(file = paste(folder.output,"/",pdfname,sep=""),paper="a4r",width=11,height=8)
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  treatmentids<-expCond$treatmentNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  nFactors<-inputparameters$nFactor
  namesFactors<-inputParameters$namesFactors
  levFactor<-inputparameters$levFactor
  prenormTable[,1:nTF]<-apply(prenormTable[,1:nTF],c(1,2),function(x) {r<-ifelse(x==0,NA,x)})
  prenormTable[,"group"]<-rep(rep(seq(1,nrepeat*ntreatment,by=1),ntime),nexp)
  # Linear scale
  if (nFactors==1){
	for(i in 1:nTF){
		prenormTabletmp<-prenormTable
		prenormTabletmp[,"TF"]<-prenormTable[,i]
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		NA0<-apply(data.frame(is.na(prenormTabletmp[,"TF"]),prenormTabletmp[,"TF"]==0),
			1,function(x){
			rr<-ifelse(any(x==TRUE),TRUE,FALSE)
			return(rr)
		})
		if(all(NA0)){
			next
		}
		summaryData<-ddply(prenormTabletmp,.(experiment,time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=prenormTabletmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			facet_wrap( ~ experiment, ncol=nCol)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab("Raw Intensity") +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else {
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(prenormTable))
		for (kk in 1:levFactor[jj]){
			prenormTabletmp<-prenormTable[
				which(prenormTable[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				prenormTabletmp[,"TF"]<-prenormTabletmp[,i]
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				titlei<-tf_ids[i]
				NA0<-apply(data.frame(is.na(prenormTabletmp[,"TF"]),prenormTabletmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				summaryData<-ddply(prenormTabletmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=prenormTabletmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab("Raw Intensity") +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
   	}
  }
  dev.off()
  # Log scale
  pdf(file = paste(folder.output,"/","Log",pdfname,sep=""),paper="a4r",width=11,height=8)
  prenormTablelog<-data.frame(
	apply(prenormTable[,1:nTF],c(1,2),function(x){
		r<-ifelse((!is.na(x)||x!=0),log(x,2),x)}),
		prenormTable[,(nTF+1):ncol(prenormTable)])
  if (nFactors==1){
	for(i in 1:nTF){
		prenormTabletmp<-prenormTablelog
		prenormTabletmp[,"TF"]<-prenormTabletmp[,i]
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		NA0<-apply(data.frame(is.na(prenormTabletmp[,"TF"]),prenormTabletmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		summaryData<-ddply(prenormTabletmp,.(experiment,time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=prenormTabletmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			facet_wrap( ~ experiment, ncol=nCol)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab("Raw Intensity (log scale)") +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else {
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(prenormTable))
		for (kk in 1:levFactor[jj]){
			prenormTabletmp<-prenormTablelog[
				which(prenormTablelog[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				prenormTabletmp[,"TF"]<-prenormTabletmp[,i]
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				titlei<-tf_ids[i]
				NA0<-apply(data.frame(is.na(prenormTabletmp[,"TF"]),prenormTabletmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				summaryData<-ddply(prenormTabletmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=prenormTabletmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab("Raw Intensity (log scale)") +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
   	}
  }
   dev.off()
}

# plotWithInArrayShiftingData
#============================================================================================

plotWithInArrayShiftingData<-function(dataAll=dataAll,inputparameters=inputParameters,
	expCond=expCond,comp=comp){
  folder.output=inputparameters$folder.output
  library(ggplot2)
  library(reshape)
  library (RColorBrewer)
    library(plyr)
  withInArray<-dataAll$dataWithInShift
  if (comp=="Gluc"){
	pdfname<-"WithInArrayShiftingDataGluc.pdf"
  } else {
	pdfname<-"WithInArrayShiftingDataFluc.pdf"
  }
  #Linear scale
  pdf(file = paste(folder.output,"/",pdfname,sep=""),paper="a4r",width=11,height=8)
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  treatmentids<-expCond$treatmentNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  nFactors<-inputparameters$nFactor
  namesFactors<-inputParameters$namesFactors
  levFactor<-inputparameters$levFactor
  withInArray[,1:nTF]<-apply(withInArray[,1:nTF],c(1,2),function(x) {r<-ifelse(x==0,NA,x)})
  withInArray[,"group"]<-rep(rep(seq(1,nrepeat*ntreatment,by=1),ntime),nexp)
  if (nFactors==1){
	for(i in 1:nTF){
		withInArraytmp<-withInArray
		withInArraytmp[,"TF"]<-withInArray[,i]
		NA0<-apply(data.frame(is.na(withInArraytmp[,"TF"]),withInArraytmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		summaryData<-ddply(withInArraytmp,.(experiment,time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=withInArraytmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			facet_wrap( ~ experiment, ncol=nCol)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab("Raw Intensity") +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else {
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(withInArray))
		for (kk in 1:levFactor[jj]){
			withInArraytmp<-withInArray[which(withInArray[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				withInArraytmp[,"TF"]<-withInArraytmp[,i]
				NA0<-apply(data.frame(is.na(withInArraytmp[,"TF"]),withInArraytmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				titlei<-tf_ids[i]
				summaryData<-ddply(withInArraytmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=withInArraytmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab("Raw Intensity") +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
   }
  }
   dev.off()
  # Logarithmic scale
  pdf(file = paste(folder.output,"/","Log",pdfname,sep=""),paper="a4r",width=11,height=8)
  withInArraylog<-data.frame(
	apply(withInArray[,1:nTF],c(1,2),function(x){
		r<-ifelse((!is.na(x)||x!=0),log(x,2),x)}),
		withInArray[,(nTF+1):ncol(withInArray)])
  if (nFactors==1){
	for(i in 1:nTF){
		withInArraytmp<-withInArraylog
		withInArraytmp[,"TF"]<-withInArraylog[,i]
		NA0<-apply(data.frame(is.na(withInArraytmp[,"TF"]),withInArraytmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		summaryData<-ddply(withInArraytmp,.(experiment,time,treatment),summarise, 
				median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
				sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
				nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=withInArraytmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			facet_wrap( ~ experiment, ncol=nCol)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab("Raw Intensity (log scale)") +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else {
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(withInArray))
		for (kk in 1:levFactor[jj]){
			withInArraytmp<-withInArraylog[which(withInArraylog[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				withInArraytmp[,"TF"]<-withInArraytmp[,i]
				NA0<-apply(data.frame(is.na(withInArraytmp[,"TF"]),withInArraytmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				titlei<-tf_ids[i]
				summaryData<-ddply(withInArraytmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=withInArraytmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab("Raw Intensity (log scale)") +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
   }
  }
   dev.off()
}

# plotnormTA
#============================================================================================

plotnormTA<-function(dataAll=dataFluc,inputparameters=inputParameters,
	expCond=expCond){
  folder.output=inputparameters$folder.output
  library(ggplot2)
  library(reshape)
  library (RColorBrewer)
    library(plyr)
  normTA<-dataAll$dataTAnormind
  pdfname<-"normTAFluc.pdf"
  pdf(file = paste(folder.output,"/",pdfname,sep=""),paper="a4r",width=11,height=8)
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  treatmentids<-expCond$treatmentNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  nFactors<-inputparameters$nFactor
  namesFactors<-inputParameters$namesFactors
  levFactor<-inputparameters$levFactor
  normTA[,1:nTF]<-apply(normTA[,1:nTF],c(1,2),function(x) {r<-ifelse(x==0,NA,x)})
  normTA[,"group"]<-rep(rep(seq(1,nrepeat*ntreatment,by=1),ntime),nexp)
  if (nFactors==1){
	for(i in 1:nTF){
		normTAtmp<-normTA
		normTAtmp[,"TF"]<-normTA[,i]
		NA0<-apply(data.frame(is.na(normTAtmp[,"TF"]),normTAtmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		summaryData<-ddply(normTAtmp,.(experiment,time,treatment),summarise, 
			median=ifelse(median(TF,na.rm=TRUE)=="NaN",NA, mean(TF,na.rm=TRUE)),
			sd=ifelse(sd(TF,na.rm=TRUE)=="NaN",NA, sd(TF,na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF))==0,1,sum(!is.na(TF))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=normTAtmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			facet_wrap( ~ experiment, ncol=nCol)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab("TA Normalized Intensity") +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else{
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(normTA))
		for (kk in 1:levFactor[jj]){
			normTAtmp<-normTA[which(normTA[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				titlei<-tf_ids[i]
				normTAtmp[,"TF"]<-normTAtmp[,i]
				NA0<-apply(data.frame(is.na(normTAtmp[,"TF"]),normTAtmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				summaryData<-ddply(normTAtmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF,na.rm=TRUE)=="NaN",NA, mean(TF,na.rm=TRUE)),
					sd=ifelse(sd(TF,na.rm=TRUE)=="NaN",NA, sd(TF,na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF))==0,1,sum(!is.na(TF))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=normTAtmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab("TA Normalized Intensity") +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
	}
   }	  
   dev.off()
   # Logarithmic scale
    pdf(file = paste(folder.output,"/","Log",pdfname,sep=""),paper="a4r",width=11,height=8)
    normTAlog<-data.frame(
		apply(normTA[,1:nTF],c(1,2),function(x){
		r<-ifelse((!is.na(x)||x!=0),log(x,2),x)}),
		normTA[,(nTF+1):ncol(normTA)])
	if (nFactors==1){
		for(i in 1:nTF){
			normTAtmp<-normTAlog
			normTAtmp[,"TF"]<-normTAlog[,i]
			NA0<-apply(data.frame(is.na(normTAtmp[,"TF"]),normTAtmp[,"TF"]==0),
				1,function(x){
					rr<-ifelse(any(x==TRUE),TRUE,FALSE)
					return(rr)
			})
			if(all(NA0)){
				next
			}
			nRow=ceiling(-1+sqrt(2+nexp))
			nCol=2+nRow
			titlei<-tf_ids[i]
			summaryData<-ddply(normTAtmp,.(experiment,time,treatment),summarise, 
				median=ifelse(median(TF,na.rm=TRUE)=="NaN",NA, mean(TF,na.rm=TRUE)),
			sd=ifelse(sd(TF,na.rm=TRUE)=="NaN",NA, sd(TF,na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF))==0,1,sum(!is.na(TF))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=normTAtmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			facet_wrap( ~ experiment, ncol=nCol)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab("TA Normalized Intensity (log scale)") +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else{
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(normTA))
		for (kk in 1:levFactor[jj]){
			normTAtmp<-normTAlog[
				which(normTAlog[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				normTAtmp[,"TF"]<-normTAtmp[,i]
				NA0<-apply(data.frame(is.na(normTAtmp[,"TF"]),normTAtmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				titlei<-tf_ids[i]
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				summaryData<-ddply(normTAtmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF,na.rm=TRUE)=="NaN",NA, mean(TF,na.rm=TRUE)),
					sd=ifelse(sd(TF,na.rm=TRUE)=="NaN",NA, sd(TF,na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF))==0,1,sum(!is.na(TF))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=normTAtmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab("TA Normalized Intensity (log scale)") +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
	}
   }	  
   dev.off()
}

# plotnormTAvsControl
#============================================================================================

plotnormTAvsControl<-function(dataAll=dataFluc,inputparameters=inputParameters,
	expCond=expCond){
  folder.output=inputparameters$folder.output
  library(ggplot2)
  library(reshape)
  library (RColorBrewer)
    library(plyr)
  normTA<-dataAll$normTAvsCtrl
  pdfname<-"normTAvsControlFluc.pdf"
  pdf(file = paste(folder.output,"/",pdfname,sep=""),paper="a4r",width=11,height=8)
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  treatmentids<-expCond$treatmentNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  nFactors<-inputparameters$nFactor
  namesFactors<-inputParameters$namesFactors
  levFactor<-inputparameters$levFactor
  normTA[,1:nTF]<-apply(normTA[,1:nTF],c(1,2),function(x) {r<-ifelse(x==0,NA,x)})
  normTA[,"group"]<-rep(rep(seq(1,nrepeat*ntreatment,by=1),ntime),nexp)
  if (nFactors==1){
	for(i in 1:nTF){
		normTAtmp<-normTA
		normTAtmp[,"TF"]<-normTA[,i]
		NA0<-apply(data.frame(is.na(normTAtmp[,"TF"]),normTAtmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		summaryData<-ddply(normTAtmp,.(experiment,time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=normTAtmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			facet_wrap( ~ experiment, ncol=nCol)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab("TA Normalized versus Control Intensity") +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else{
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(normTA))
		for (kk in 1:levFactor[jj]){
			normTAtmp<-normTA[which(normTA[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				normTAtmp[,"TF"]<-normTAtmp[,i]
				NA0<-apply(data.frame(is.na(normTAtmp[,"TF"]),normTAtmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				summaryData<-ddply(normTAtmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=normTAtmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab("TA Normalized versus Control Intensity") +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
	}
   }	  
   dev.off()
   # Logarithmic scale
    pdf(file = paste(folder.output,"/","Log",pdfname,sep=""),paper="a4r",width=11,height=8)
    normTAlog<-data.frame(
		apply(normTA[,1:nTF],c(1,2),function(x){
		r<-ifelse((!is.na(x)||x!=0),log(x,2),x)}),
		normTA[,(nTF+1):ncol(normTA)])
	if (nFactors==1){
		for(i in 1:nTF){
			normTAtmp<-normTAlog
			normTAtmp[,"TF"]<-normTAlog[,i]
			NA0<-apply(data.frame(is.na(normTAtmp[,"TF"]),normTAtmp[,"TF"]==0),
				1,function(x){
					rr<-ifelse(any(x==TRUE),TRUE,FALSE)
					return(rr)
			})
			if(all(NA0)){
				next
			}
			nRow=ceiling(-1+sqrt(2+nexp))
			nCol=2+nRow
			titlei<-tf_ids[i]
			summaryData<-ddply(normTAtmp,.(experiment,time,treatment),summarise, 
				median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
				sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
				nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
			summaryData$TF<-summaryData$median
			summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
			summaryData$lcl<-summaryData$TF-1.96*summaryData$se
			summaryData$ucl<-summaryData$TF+1.96*summaryData$se
			if (all(!is.finite(summaryData$TF))){
				next
			}
			pp<-ggplot(data=normTAtmp,
				aes(x= time,y= TF,colour= treatment))+
				geom_point(na.rm=TRUE)+
				facet_wrap( ~ experiment, ncol=nCol)+
				geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
				scale_fill_brewer(type="qual",palette=1,
					name = "Treatment")+
				xlab("Time") +
				ylab("TA Normalized versus Control Intensity (log scale)") +
				labs(title =titlei)+
				scale_colour_brewer(type="qual",palette=1,
					name = "Treatment")+
				theme(panel.background = element_rect(fill="transparent"))+
				theme(strip.background=element_rect(fill="white"))+
				theme(axis.ticks=element_line(colour="black"))+
				theme(axis.title=element_text(size=16))+
				theme(panel.border=element_rect(colour="black",fill="transparent"))
			print(pp)
		}
  } else{
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(normTA))
		for (kk in 1:levFactor[jj]){
			normTAtmp<-normTAlog[
				which(normTAlog[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				normTAtmp[,"TF"]<-normTAtmp[,i]
				NA0<-apply(data.frame(is.na(normTAtmp[,"TF"]),normTAtmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				summaryData<-ddply(normTAtmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=normTAtmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab("TA Normalized versus Control Intensity(log scale)") +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
	}
   }	  
   dev.off()
}

# plotNormalization
#============================================================================================
 
plotNormalization<-function(dataAll=dataAllFluc,inputparameters=inputParameters,expCond=expCond){
	folder.output=inputparameters$folder.output
	library(ggplot2)
	library(reshape)
	library (RColorBrewer)
	library(Hmisc)
	library(plyr)
	normdata<-dataAll[["dataScaleFree"]]
	ntreatment<-inputparameters$ntreat
	tf_ids<-expCond$tfNames
	ntime<-inputparameters$ntime
	timevector<-inputparameters$timeVector
	treatmentids<-expCond$treatmentNames
	nTF<-inputparameters$nTF
	nexp<-inputparameters$nexp
	nrepeat<-max(normdata[,"repeat"],na.rm=TRUE)
	nFactors<-inputparameters$nFactor
	namesFactors<-inputParameters$namesFactors
	levFactor<-inputparameters$levFactor
	plotnorm(dataAll=normdata,inputparameters=inputParameters,
		expCond=expCond,pdfname="Normalization.pdf",ytitle="Normalized Intensity")
	plotnormAverage(dataAll=normdata ,inputparameters=inputParameters,
		expCond=expCond,pdfname="NormalizationAverage.pdf",ytitle="Average Normalized Intensity")
	if (nTF<25){
		plotnormAverageSmallArrays(dataAll=normdata ,inputparameters=inputParameters,
		expCond=expCond,pdfname="NormalizationAverageAllTFtogether.pdf",
		ytitle="Average Normalized Intensity")
	}
	if (any(inputparameters$nameUntreated!="none")){
		plotnorm(dataAll=dataAll$normalizedvsCtrl,
			inputparameters=inputParameters,expCond=expCond,
			pdfname="NormalizationvsControl.pdf",
			ytitle="Normalized Intensity versus Control")
		plotnormAverage(dataAll=dataAll$normalizedvsCtrl,
			inputparameters=inputParameters,expCond=expCond,
			pdfname="NormalizationAveragevsControl.pdf",
			ytitle="Average Normalized Intensity versus Control")
		if (nTF<25){
			plotnormAverageSmallArrays(dataAll=dataAll$normalizedvsCtrl,
				inputparameters=inputParameters,expCond=expCond,
				pdfname="NormalizationAverageAllTFtogethervsControl.pdf",
				ytitle="Average Normalized Intensity versus Control")
		}
	}
}
	
# plotnorm
#============================================================================================

plotnorm<-function(dataAll=datat0Norm,inputparameters=inputParameters,
	expCond=expCond,pdfname=pdfname,ytitle="Normalized Intensity"){
  folder.output=inputparameters$folder.output	
  library(ggplot2)
  library(reshape)
  library (RColorBrewer)
    library(plyr)
  norm<-dataAll
  #Linear scale
  pdf(file = paste(folder.output,"/",pdfname,sep=""),paper="a4r",width=11,height=8)
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  treatmentids<-expCond$treatmentNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  nFactors<-inputparameters$nFactor
  namesFactors<-inputParameters$namesFactors
  levFactor<-inputparameters$levFactor
  norm[,1:nTF]<-apply(norm[,1:nTF],c(1,2),function(x) {r<-ifelse(x==0,NA,x)})
  norm[,"group"]<-rep(rep(seq(1,nrepeat*ntreatment,by=1),ntime),nexp)
  if (nFactors==1){
	for(i in 1:nTF){
		normtmp<-norm
		normtmp[,"TF"]<-norm[,i]
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		summaryData<-ddply(normtmp,.(experiment,time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=normtmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			facet_wrap( ~ experiment, ncol=nCol)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab(ytitle) +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
			name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else {
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(norm))
		for (kk in 1:levFactor[jj]){
			normtmp<-norm[which(norm[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				normtmp[,"TF"]<-normtmp[,i]
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				titlei<-tf_ids[i]
				NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				summaryData<-ddply(normtmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=normtmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab(ytitle) +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
	}
   }	
   dev.off()
   # Log scale
   pdf(file = paste(folder.output,"/","Log",pdfname,sep=""),paper="a4r",width=11,height=8)
   normlog<-data.frame(
		apply(norm[,1:nTF],c(1,2),function(x){
		r<-ifelse((!is.na(x)||x!=0),log(x,2),x)}),
		norm[,(nTF+1):ncol(norm)])
   if (nFactors==1){
	for(i in 1:nTF){
		normtmp<-normlog
		normtmp[,"TF"]<-normlog[,i]
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		summaryData<-ddply(normtmp,.(experiment,time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=normtmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			facet_wrap( ~ experiment, ncol=nCol)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab(paste(ytitle," (log scale)",sep="")) +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else {
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(norm))
		for (kk in 1:levFactor[jj]){
			normtmp<-normlog[
				which(normlog[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				normtmp[,"TF"]<-normtmp[,i]
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				titlei<-tf_ids[i]
				NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				summaryData<-ddply(normtmp,.(experiment,time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=normtmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					facet_wrap( ~ experiment, ncol=nCol)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab(paste(ytitle," (log scale)",sep="")) +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
	}
   }	
   dev.off()
}

# plotnormAverage
#============================================================================================

plotnormAverage<-function(dataAll=datat0Norm,inputparameters=inputParameters,
	expCond=expCond,pdfname=pdfname,ytitle="Normalized Intensity (average)"){
  folder.output=inputparameters$folder.output	
  library(ggplot2)
  library(reshape)
  library (RColorBrewer)
  library(plyr)
  norm<-dataAll
  #Linear scale
  pdf(file = paste(folder.output,"/",pdfname,sep=""),paper="a4r",width=11,height=8)
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  treatmentids<-expCond$treatmentNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  nFactors<-inputparameters$nFactor
  namesFactors<-inputParameters$namesFactors
  levFactor<-inputparameters$levFactor
  norm[,1:nTF]<-apply(norm[,1:nTF],c(1,2),function(x) {r<-ifelse(x==0,NA,x)})
  norm[,"group"]<-rep(rep(seq(1,nrepeat*ntreatment,by=1),ntime),nexp)
  if (nFactors==1){
	for(i in 1:nTF){
		normtmp<-norm
		normtmp[,"TF"]<-norm[,i]
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		summaryData<-ddply(normtmp,.(time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=normtmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab(ytitle) +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
			name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else {
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(norm))
		for (kk in 1:levFactor[jj]){
			normtmp<-norm[which(norm[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				normtmp[,"TF"]<-normtmp[,i]
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				titlei<-tf_ids[i]
				NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
					1,function(x){
					rr<-ifelse(any(x==TRUE),TRUE,FALSE)
					return(rr)
				})
				if(all(NA0)){
					next
				}
				summaryData<-ddply(normtmp,.(time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=normtmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab(ytitle) +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
	}
   }	
   dev.off()
   # Log scale
   pdf(file = paste(folder.output,"/","Log",pdfname,sep=""),paper="a4r",width=11,height=8)
   normlog<-data.frame(
		apply(norm[,1:nTF],c(1,2),function(x){
		r<-ifelse((!is.na(x)||x!=0),log(x,2),x)}),
		norm[,(nTF+1):ncol(norm)])
   if (nFactors==1){
	for(i in 1:nTF){
		normtmp<-normlog
		normtmp[,"TF"]<-normlog[,i]
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		summaryData<-ddply(normtmp,.(time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
		summaryData$TF<-summaryData$median
		summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
		summaryData$lcl<-summaryData$TF-1.96*summaryData$se
		summaryData$ucl<-summaryData$TF+1.96*summaryData$se
		if (all(!is.finite(summaryData$TF))){
			next
		}
		pp<-ggplot(data=normtmp,
			aes(x= time,y= TF,colour= treatment))+
			geom_point(na.rm=TRUE)+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			scale_fill_brewer(type="qual",palette=1,
				name = "Treatment")+
			xlab("Time") +
			ylab(paste(ytitle," (log scale)",sep="")) +
			labs(title =titlei)+
			scale_colour_brewer(type="qual",palette=1,
				name = "Treatment")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=16))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		print(pp)
	}
  } else {
	for (jj in 1:nFactors){
		posjj<-grep(namesFactors[jj],colnames(norm))
		for (kk in 1:levFactor[jj]){
			normtmp<-normlog[
				which(normlog[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
			for(i in 1:nTF){
				normtmp[,"TF"]<-normtmp[,i]
				nRow=ceiling(-1+sqrt(2+nexp))
				nCol=2+nRow
				titlei<-tf_ids[i]
				NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
					1,function(x){
						rr<-ifelse(any(x==TRUE),TRUE,FALSE)
						return(rr)
				})
				if(all(NA0)){
					next
				}
				summaryData<-ddply(normtmp,.(time,treatment),summarise, 
					median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
					sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
					nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
				summaryData$TF<-summaryData$median
				summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
				summaryData$lcl<-summaryData$TF-1.96*summaryData$se
				summaryData$ucl<-summaryData$TF+1.96*summaryData$se
				if (all(!is.finite(summaryData$TF))){
					next
				}
				pp<-ggplot(data=normtmp,
					aes(x= time,y= TF,colour= treatment))+
					geom_point(na.rm=TRUE)+
					geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
					scale_fill_brewer(type="qual",palette=1,
						name = "Treatment")+
					xlab("Time") +
					ylab(paste(ytitle," (log scale)",sep="")) +
					labs(title =titlei)+
					scale_colour_brewer(type="qual",palette=1,
						name = "Treatment")+
					theme(panel.background = element_rect(fill="transparent"))+
					theme(strip.background=element_rect(fill="white"))+
					theme(axis.ticks=element_line(colour="black"))+
					theme(axis.title=element_text(size=16))+
					theme(panel.border=element_rect(colour="black",fill="transparent"))
				print(pp)
			}
		}
	}
   }	
   dev.off()
}
# plotnormAverageSmallArrays
#============================================================================================

plotnormAverageSmallArrays<-function(dataAll=datat0Norm,inputparameters=inputParameters,
	expCond=expCond,pdfname=pdfname,ytitle="Normalized Intensity (average)"){
   pdf(file = paste(folder.output,"/",pdfname,sep=""),paper="a4r",width=11,height=8)
   ntreatment<-inputparameters$ntreat
     library(plyr)
   tf_ids<-expCond$tfNames
   tfBlank<-inputparameters$tfBlank
   nTF<-inputparameters$nTF
   normdata<-dataAll
   if (!is.null(tfBlank)){
    if (tfBlank%in%tf_ids){
		tf_ids<-tf_ids[-which(tf_ids==tfBlank)]
		nTF=nTF-1
		normdata<-normdata[,-grep(tfBlank,colnames(normdata))]
	}
   }
   ntime<-inputparameters$ntime
   timevector<-inputparameters$timeVector
   treatmentids<-expCond$treatmentNames
   nexp<-inputparameters$nexp
   nrepeat<-inputparameters$nrepeat
   nFactors<-inputparameters$nFactor
   namesFactors<-inputParameters$namesFactors
   levFactor<-inputparameters$levFactor
   if (nFactors==1){
		n = sqrt(nTF+1)
		nCol = ceiling(n)
		nRow = ceiling(n)
		library(grid)
		vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
		normtmp<-normdata
		for (ii in 1:nTF){
			normtmp[,"TF"]<-normdata[,ii]
			titlei<-tf_ids[ii]
			library(plyr)
			NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
				1,function(x){
					rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
			})
			if(all(NA0)){
				next
			}
			if (abs(median(normtmp[,"TF"],na.rm=TRUE))<1){
				normtmp[,"TF"]*100
				ytitle="Percentage Normalized Intensity (average)"
			}
			summaryData<-ddply(normtmp,.(time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
			summaryData$TF<-summaryData$median
			summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
			summaryData$lcl<-summaryData$TF-1.96*summaryData$se
			summaryData$ucl<-summaryData$TF+1.96*summaryData$se
			if (all(!is.finite(summaryData$TF))){
				next
			}
			ymin<-min(normtmp[,"TF"],na.rm=TRUE)*0.9
			ymax<-max(normtmp[,"TF"],na.rm=TRUE)*0.9
			ymean<-0.5*(ymin+ymax)
			dodge <- position_dodge(width=0.9)
			pp<-ggplot(data=normtmp,
				aes(x= time,y= TF,colour= treatment))+
				geom_point(na.rm=TRUE)+
				geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
				scale_fill_brewer(type="qual",palette=1,
					name = "Treatment")+
				xlab("")+
				ylab("")+
				labs(title=titlei)+
				scale_colour_brewer(type="qual",palette=1,
					name = "Treatment")+
				theme(panel.background = element_rect(fill="transparent"))+
				theme(strip.background=element_rect(fill="white"))+
				theme(axis.ticks=element_line(colour="black"))+
				theme(axis.title=element_text(size=12))+
				theme(panel.border=element_rect(colour="black",fill="transparent"))+
				theme(strip.text=element_text(size=12,face="bold"))+
				theme(axis.text=element_text(size=12,colour="black",face="bold"))+
				theme(legend.position="none")+
				theme(plot.title=element_text(size=12, colour="black",face="bold"))+
				scale_y_continuous(breaks=c(ymin,ymean,ymax),
					labels=c(round(ymin,1),round(ymean,1),round(ymax,1)))
				assign(paste("plot_",tf_ids[ii],sep=""),pp)
		}
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nRow, nCol)))
		idtf<-0
		for (ll in 1:nRow){
			for (mm in 1:nCol){
				idtf<-idtf+1
				if (idtf<=nTF){
					try(print(get(paste("plot_",tf_ids[idtf],sep="")), 
						vp = vplayout(ll, mm)),TRUE)
				}
			}			
		}
	} else {
		for (jj in 1:nFactors){
			posjj<-grep(namesFactors[jj],colnames(normdata))
			for (kk in 1:levFactor[jj]){
				n = sqrt(nTF+1)
				nCol = ceiling(n)
				nRow = ceiling(n)
				library(grid)
				vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
				normtmp<-normdata[which(normdata[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
				for (ii in 1:nTF){
					normtmp[,"TF"]<-normdata[,ii]
					titlei<-tf_ids[ii]
					library(plyr)
					NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
						1,function(x){
							rr<-ifelse(any(x==TRUE),TRUE,FALSE)
							return(rr)
					})
					if(all(NA0)){
						next
					}
					if (abs(median(normtmp[,"TF"],na.rm=TRUE))<1){
						normtmp[,"TF"]*100
						ytitle="Percentage Normalized Intensity (average)"
					}
					summaryData<-ddply(normtmp,.(time,treatment),summarise, 
						median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
						sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
						nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
					summaryData$TF<-summaryData$median
					summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
					summaryData$lcl<-summaryData$TF-1.96*summaryData$se
					summaryData$ucl<-summaryData$TF+1.96*summaryData$se
					if (all(!is.finite(summaryData$TF))){
						next
					}
					ymin<-min(normtmp[,"TF"],na.rm=TRUE)*0.9
					ymax<-max(normtmp[,"TF"],na.rm=TRUE)*0.9
					ymean<-0.5*(ymin+ymax)
					dodge <- position_dodge(width=0.9)
					pp<-ggplot(data=normtmp,
						aes(x= time,y= TF,colour= treatment))+
						geom_point(na.rm=TRUE)+
						geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
						scale_fill_brewer(type="qual",palette=1,
							name = "Treatment")+
						xlab("")+
						ylab("")+
						labs(title=titlei)+
						scale_colour_brewer(type="qual",palette=1,
							name = "Treatment")+
						theme(panel.background = element_rect(fill="transparent"))+
						theme(strip.background=element_rect(fill="white"))+
						theme(axis.ticks=element_line(colour="black"))+
						theme(axis.title=element_text(size=12))+
						theme(panel.border=element_rect(colour="black",fill="transparent"))+
						theme(strip.text=element_text(size=12,face="bold"))+
						theme(axis.text=element_text(size=12,colour="black",face="bold"))+
						theme(legend.position="none")+
						theme(plot.title=element_text(size=12, colour="black",face="bold"))+
						scale_y_continuous(breaks=c(ymin,ymean,ymax),
							labels=c(round(ymin,1),round(ymean,1),round(ymax,1)))
						assign(paste("plot_",tf_ids[ii],sep=""),pp)
				}
				grid.newpage()
				pushViewport(viewport(layout = grid.layout(nRow, nCol)))
				idtf<-0
				for (ll in 1:nRow){
					for (mm in 1:nCol){
						idtf<-idtf+1
						if (idtf<=nTF){
							try(print(get(paste("plot_",tf_ids[idtf],sep="")), 
								vp = vplayout(ll, mm)),TRUE)
						}
					}			
				}
			}
		}
	}
	dev.off()	
	# In logartimic scale
	pdf(file = paste(folder.output,"/Log",pdfname,sep=""),paper="a4r",width=11,height=8)
	if (nFactors==1){
		n = sqrt(nTF+1)
		nCol = ceiling(n)
		nRow = ceiling(n)
		library(grid)
		vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
		normdatalog<-apply(normdata[,1:nTF],c(1,2),function(x){
			x<-as.numeric(x)
			if (!is.na(x) && x!=0){
				rr<-log(x,2)
			} else {
				rr<-ifelse(is.na(x),NA,0)
			}
			return(rr)
		})	
		normdatalog<-data.frame(normdatalog,normdata[,((nTF+1):ncol(normdata))])
		normtmp<-normdatalog
		for (ii in 1:nTF){
			normtmp[,"TF"]<-normdatalog[,ii]
			titlei<-tf_ids[ii]
			library(plyr)
			NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
				1,function(x){
					rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
			})
			if(all(NA0)){
				next
			}
			if (abs(median(normtmp[,"TF"],na.rm=TRUE))<1){
				normtmp[,"TF"]*100
				ytitle="Percentage Log2 Normalized Intensity (average)"
			}
			summaryData<-ddply(normtmp,.(time,treatment),summarise, 
			median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
			sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
			nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
			summaryData$TF<-summaryData$median
			summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
			summaryData$lcl<-summaryData$TF-1.96*summaryData$se
			summaryData$ucl<-summaryData$TF+1.96*summaryData$se
			if (all(!is.finite(summaryData$TF))){
				next
			}
			ymin<-min(normtmp[,"TF"],na.rm=TRUE)*0.9
			ymax<-max(normtmp[,"TF"],na.rm=TRUE)*0.9
			ymean<-0.5*(ymin+ymax)
			dodge <- position_dodge(width=0.9)
			pp<-ggplot(data=normtmp,
				aes(x= time,y= TF,colour= treatment))+
				geom_point(na.rm=TRUE)+
				geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
				scale_fill_brewer(type="qual",palette=1,
					name = "Treatment")+
				xlab("")+
				ylab("")+
				labs(title=titlei)+
				scale_colour_brewer(type="qual",palette=1,
					name = "Treatment")+
				theme(panel.background = element_rect(fill="transparent"))+
				theme(strip.background=element_rect(fill="white"))+
				theme(axis.ticks=element_line(colour="black"))+
				theme(axis.title=element_text(size=12))+
				theme(panel.border=element_rect(colour="black",fill="transparent"))+
				theme(strip.text=element_text(size=12,face="bold"))+
				theme(axis.text=element_text(size=12,colour="black",face="bold"))+
				theme(legend.position="none")+
				theme(plot.title=element_text(size=12, colour="black",face="bold"))+
				scale_y_continuous(breaks=c(ymin,ymean,ymax),
					labels=c(round(ymin,1),round(ymean,1),round(ymax,1)))
				assign(paste("plot_",tf_ids[ii],sep=""),pp)
		}
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nRow, nCol)))
		idtf<-0
		for (ll in 1:nRow){
			for (mm in 1:nCol){
				idtf<-idtf+1
				if (idtf<=nTF){
					try(print(get(paste("plot_",tf_ids[idtf],sep="")), 
						vp = vplayout(ll, mm)),TRUE)
				}
			}			
		}
	} else {
		for (jj in 1:nFactors){
			posjj<-grep(namesFactors[jj],colnames(normdata))
			for (kk in 1:levFactor[jj]){
				n = sqrt(nTF+1)
				nCol = ceiling(n)
				nRow = ceiling(n)
				library(grid)
				vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
				normtmp<-normdatalog[which(normdatalog[,posjj]==expCond$levelsNames[[namesFactors[jj]]][kk]),]
				for (ii in 1:nTF){
					normtmp[,"TF"]<-normdatalog[,ii]
					titlei<-tf_ids[ii]
					library(plyr)
					NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
						1,function(x){
							rr<-ifelse(any(x==TRUE),TRUE,FALSE)
							return(rr)
					})
					if(all(NA0)){
						next
					}
					if (abs(median(normtmp[,"TF"],na.rm=TRUE))<1){
						normtmp[,"TF"]*100
						ytitle="Percentage Log2 Normalized Intensity (average)"
					}
					summaryData<-ddply(normtmp,.(time,treatment),summarise, 
						median=ifelse(median(TF[TF!=0],na.rm=TRUE)=="NaN",NA, mean(TF[TF!=0],na.rm=TRUE)),
						sd=ifelse(sd(TF[TF!=0],na.rm=TRUE)=="NaN",NA, sd(TF[TF!=0],na.rm=TRUE)),
						nrep=ifelse(sum(!is.na(TF[TF!=0]))==0,1,sum(!is.na(TF[TF!=0]))))
					summaryData$TF<-summaryData$median
					summaryData$se<-summaryData$sd/(summaryData$nrep)^0.5
					summaryData$lcl<-summaryData$TF-1.96*summaryData$se
					summaryData$ucl<-summaryData$TF+1.96*summaryData$se
					if (all(!is.finite(summaryData$TF))){
						next
					}
					ymin<-min(normtmp[,"TF"],na.rm=TRUE)*0.9
					ymax<-max(normtmp[,"TF"],na.rm=TRUE)*0.9
					ymean<-0.5*(ymin+ymax)
					dodge <- position_dodge(width=0.9)
					pp<-ggplot(data=normtmp,
						aes(x= time,y= TF,colour= treatment))+
						geom_point(na.rm=TRUE)+
						geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
						scale_fill_brewer(type="qual",palette=1,
							name = "Treatment")+
						xlab("")+
						ylab("")+
						labs(title=titlei)+
						scale_colour_brewer(type="qual",palette=1,
							name = "Treatment")+
						theme(panel.background = element_rect(fill="transparent"))+
						theme(strip.background=element_rect(fill="white"))+
						theme(axis.ticks=element_line(colour="black"))+
						theme(axis.title=element_text(size=12))+
						theme(panel.border=element_rect(colour="black",fill="transparent"))+
						theme(strip.text=element_text(size=12,face="bold"))+
						theme(axis.text=element_text(size=12,colour="black",face="bold"))+
						theme(legend.position="none")+
						theme(plot.title=element_text(size=12, colour="black",face="bold"))+
						scale_y_continuous(breaks=c(ymin,ymean,ymax),
							labels=c(round(ymin,1),round(ymean,1),round(ymax,1)))
						assign(paste("plot_",tf_ids[ii],sep=""),pp)
				}
				grid.newpage()
				pushViewport(viewport(layout = grid.layout(nRow, nCol)))
				idtf<-0
				for (ll in 1:nRow){
					for (mm in 1:nCol){
						idtf<-idtf+1
						if (idtf<=nTF){
							try(print(get(paste("plot_",tf_ids[idtf],sep="")), 
								vp = vplayout(ll, mm)),TRUE)
						}
					}			
				}
			}
		}
	}
	dev.off()	
}

# plotnormWeigthedAverage
#============================================================================================
 plotnormWeigthedAverage<-function(dataAll=table6ind ,inputparameters=inputparameters,
			expCond=expCond,pdfname="FinalTFAverage.pdf",ytitle="Average Normalized Intensity",
			pmatrix=pmatrix,pmatrixZ=pmatrixZ,namespmatrix=namespmatrix){ 
  folder.output=inputparameters$folder.output	
  library(ggplot2)
  library(reshape)
  library (RColorBrewer)
  library(plyr)
  norm<-dataAll
  valueColors<-inputparameters$valueColors
  #Linear scale
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  tf_ids<-gsub("/","_",tf_ids,fixed=TRUE)
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  treatmentids<-expCond$treatmentNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  nFactors<-inputparameters$nFactor
  namesFactors<-inputParameters$namesFactors
  levFactor<-inputparameters$levFactor
  norm[,1:nTF]<-apply(norm[,1:nTF],c(1,2),function(x) {r<-ifelse(x==0,NA,x)})
  norm[,"group"]<-rep(rep(seq(1,nrepeat*ntreatment,by=1),ntime),nexp)
 #tmpmatrix<-pmatrix
  #tmpmatrix[which(abs(tmpmatrix)>0.05 | is.na(tmpmatrix))]<-0.5
  #pmatrix<-tmpmatrix
  toadd<-max(-round(qnorm(min(abs(pmatrix),na.rm=TRUE)),0),4)
  if (nFactors==1){
	for(i in 1:nTF){
		wZ<-c()
		for (oo in 1:length(timevector)){
			pos_oo<-grep(paste("_T",timevector[oo],"-",sep=""),namespmatrix)
			pmatrix_oo<-pmatrix[i,pos_oo,] #identify all the p-values for a given time for all the experiments
			pmatrix_oo_Z<-pmatrixZ[i,pos_oo] #identify the z-score provided to each treatment time
			pmatrix_oo<-pmatrix_oo*matrix(sign(pmatrix_oo_Z),ncol=ncol(pmatrix_oo),
				nrow=nrow(pmatrix_oo))
			pmatrix_oo<-apply(pmatrix_oo,c(1,2),function(x){
				rr<-ifelse(abs(x)==1,0.999*sign(x),x)
				return(rr)
			})
			if (length(pos_oo)>1){
				pmatrixZ_oo<-apply(pmatrix_oo,c(1,2),function(x){
					rr<-ifelse(is.finite(x),
						ifelse(sign(x)==1,(-qnorm(x)+toadd)^2, (qnorm(abs(x))+toadd)),
						0)
					return(rr)
				})
				wZtmp<-colSums(pmatrixZ_oo)
				wZtmp<-rep(wZtmp,each=(nrepeat*length(treatmentids)))
				wZ<-c(wZ,wZtmp)
			} else{
				pmatrixZ_oo<-pmatrix_oo
				for (mm in 1:length(pmatrix_oo)){
					pmatrixZ_oo[mm]<-ifelse(is.finite(pmatrix_oo[mm]),
						ifelse(sign(pmatrix_oo[mm])==1,(-qnorm(pmatrix_oo[mm])+toadd)^2, 
							qnorm(abs(pmatrix_oo[mm]))+toadd),
						0)
				}
				wZtmp<-pmatrixZ_oo
				wZtmp<-rep(wZtmp,each=(nrepeat*length(treatmentids)))
				wZ<-c(wZ,wZtmp)
			}
		}
		orderwZ<-rep(rep(seq(1:nexp),each=(nrepeat*length(treatmentids))),ntime)
		wZ<-wZ[order(orderwZ,decreasing=FALSE)]
		normtmp<-norm
		normtmp[,"TF"]<-norm[,i]
		normtmp[,"wZ"]<-wZ
		nRow=ceiling(-1+sqrt(2+nexp))
		nCol=2+nRow
		titlei<-tf_ids[i]
		NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
			1,function(x){
				rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
		})
		if(all(NA0)){
			next
		}
		summaryData<-matrix(NA,ncol=6,nrow=(ntime*length(treatmentids)))
		idrow=0
		for (aa in 1:ntime){
			for (bb in 1:length(treatmentids)){
				idrow=idrow+1
				tmp<-normtmp[which(normtmp[,"time"]==timevector[aa] & normtmp[,"treatment"]==treatmentids[bb]),]
				if(median(tmp[which(tmp[,"TF"]!=0),"TF"],na.rm=TRUE)=="NaN"){
					mean_w=NA
				} else {
					TF<-tmp[!is.na(tmp[,"TF"]),"TF"]
					wZ<-tmp[!is.na(tmp[,"TF"]),"wZ"]
					zz<-TF*wZ
					mean_w<-ifelse(sum(wZ)!=0,sum(zz)/sum(wZ),NA)
				}
				summaryData[idrow,1]=mean_w
				nrep=ifelse(sum(!is.na(tmp[which(tmp[,"TF"]!=0),"TF"]))==0,1,
					sum(!is.na(tmp[which(tmp[,"TF"]!=0),"TF"])))
				summaryData[idrow,4]=nrep
				if(sd(tmp[which(tmp[,"TF"]!=0),"TF"],na.rm=TRUE)=="NaN"){
					sd_w=NA
				} else {
					TF<-tmp[!is.na(tmp[,"TF"]),"TF"]
					wZ<-tmp[!is.na(tmp[,"TF"]),"wZ"]
					zz<-TF*wZ
					mean_w<-ifelse(sum(wZ)!=0,sum(zz)/sum(wZ),NA)
					# Obtain the weighted variability
					nonzerowZ<-sum(which(wZ>0))
					var_w<-ifelse(sum(wZ)!=0,
						(sum(wZ*((TF-mean_w)^2))/((nonzerowZ-1)*sum(wZ)/nonzerowZ)),NA)
					sd_w=var_w^0.5
				}
				summaryData[idrow,2]=sd_w
				summaryData[idrow,3]=sd_w/nonzerowZ^0.5
				summaryData[idrow,5]=mean_w-1.96*(sd_w/nonzerowZ^0.5)
				summaryData[idrow,6]=mean_w+1.96*(sd_w/nonzerowZ^0.5)
			}
		}
		summaryData<-data.frame(time=rep(timevector,each=length(treatmentids)),treatment=rep(treatmentids,ntime),summaryData)
		colnames(summaryData)<-c("time","treatment","TF","sd_w","se","nrep","lcl","ucl")
		if (i==1){
			summaryDataWeighted<-summaryData
		} else {
			summaryDataWeighted<-rbind(summaryDataWeighted,summaryData)
		}
		if (all(!is.finite(summaryData[,"TF"]))){
			next
		}
		pp<-ggplot(data=summaryData,
			aes(x= time,y= TF,group=treatment,colour=treatment))+
			geom_line(na.rm=TRUE,colour="black")+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,fill=treatment), stat="identity")+
			xlab("Time") +
			ylab(ytitle) +
			labs(title =titlei)+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=24))+
			theme(axis.text=element_text(size=20))+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		if(is.null(valueColors)){
			pp<-pp+scale_colour_brewer(type="qual",palette=1,name = "Treatment")
				
		} else {
			pp<-pp+scale_fill_manual(name = "Treatment", values = valueColors)+
				scale_colour_manual(name = "Treatment", values = valueColors)
		}
		pdf(file = paste(folder.output,"/",tf_ids[i],"_",pdfname,sep=""),paper="a4r",width=11,height=8)
		print(pp)
		dev.off()
		# For publication
		pp<-ggplot(data=summaryData,
			aes(x= time,y= TF,group=treatment,colour=treatment))+
			geom_line(na.rm=TRUE,colour="black")+
			geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,
				fill=treatment), stat="identity")+
			xlab("") +
			ylab("") +
			labs("")+
			theme(panel.background = element_rect(fill="transparent"))+
			theme(strip.background=element_rect(fill="white"))+
			theme(axis.ticks=element_line(colour="black"))+
			theme(axis.title=element_text(size=24))+
			theme(axis.text=element_text(size=24,colour="black",face="bold"))+
			theme(legend.position="none")+
			theme(panel.border=element_rect(colour="black",fill="transparent"))
		if(is.null(valueColors)){
			pp<-pp+scale_colour_brewer(type="qual",palette=1,name = "Treatment")
				
		} else {
			pp<-pp+scale_fill_manual(name = "Treatment", values = valueColors)+
				scale_colour_manual(name = "Treatment", values = valueColors)
		}
		pdf(file = paste(folder.output,"/","forpublication_",tf_ids[i],"_",pdfname,sep=""),paper="a4r",width=11,height=8)
		print(pp)
		dev.off()
	}
  }
  summaryDataWeighted<-cbind(rep(tf_ids,each=(ntime*length(treatmentids))),summaryDataWeighted)
  write.table(summaryDataWeighted, file=paste(folder.output,"/summaryDataWeighted.txt",sep=""), sep="\t",
	row.names = FALSE,col.names = TRUE)
}

# plotnormWeightedAverageSmallArrays
#============================================================================================
plotnormWeigthedAverageSmallArrays<-function(dataAll=table6ind ,inputparameters=inputparameters,
			expCond=expCond,pdfname="FinalTFAverageAllTFtogether.pdf",ytitle="Average Normalized Intensity",
			pmatrix=pmatrix,pmatrixZ=pmatrixZ,namespmatrix=namespmatrix){
   ntreatment<-inputparameters$ntreat
     library(plyr)
   tf_ids<-expCond$tfNames
   tf_ids<-gsub("/","_",tf_ids,fixed=TRUE)
   tfBlank<-inputparameters$tfBlank
   nTF<-inputparameters$nTF
   normdata<-dataAll
   if (!is.null(tfBlank)){
    if (tfBlank%in%tf_ids){
		tf_ids<-tf_ids[-which(tf_ids==tfBlank)]
		nTF=nTF-1
		normdata<-normdata[,-grep(tfBlank,colnames(normdata))]
	}
   }
   ntime<-inputparameters$ntime
   timevector<-inputparameters$timeVector
   treatmentids<-expCond$treatmentNames
   nexp<-inputparameters$nexp
   nrepeat<-inputparameters$nrepeat
   nFactors<-inputparameters$nFactor
   namesFactors<-inputParameters$namesFactors
   levFactor<-inputparameters$levFactor
    #tmpmatrix<-pmatrix
   #tmpmatrix[which(abs(tmpmatrix)>0.05 | is.na(tmpmatrix))]<-0.5
   #pmatrix<-tmpmatrix
   toadd<-max(-round(qnorm(min(abs(pmatrix),na.rm=TRUE)),0),4)
   if (nFactors==1){
		library(grid)
		vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
		normtmp<-normdata
		for (ii in 1:nTF){
			wZ<-c()
			for (oo in 1:length(timevector)){
				pos_oo<-grep(paste("_T",timevector[oo],"-",sep=""),namespmatrix)
				pmatrix_oo<-pmatrix[ii,pos_oo,] #identify all the p-values for a given time for all the experiments
				pmatrix_oo_Z<-pmatrixZ[ii,pos_oo] #identify the z-score provided to each treatment time
				pmatrix_oo<-pmatrix_oo*matrix(sign(pmatrix_oo_Z),ncol=ncol(pmatrix_oo),
					nrow=nrow(pmatrix_oo))
				pmatrix_oo<-apply(pmatrix_oo,c(1,2),function(x){
					rr<-ifelse(abs(x)==1,0.999*sign(x),x)
					return(rr)
				})
				if (length(pos_oo)>1){
					pmatrixZ_oo<-apply(pmatrix_oo,c(1,2),function(x){
						rr<-ifelse(is.finite(x),
							ifelse(sign(x)==1,(-qnorm(x)+toadd)^2, qnorm(abs(x))+toadd),
							0)
						return(rr)
					})
					wZtmp<-colSums(pmatrixZ_oo)
					wZtmp<-rep(wZtmp,each=(nrepeat*length(treatmentids)))
					wZ<-c(wZ,wZtmp)
				} else{
					pmatrixZ_oo<-pmatrix_oo
					for (mm in 1:length(pmatrix_oo)){
						pmatrixZ_oo[mm]<-ifelse(is.finite(pmatrix_oo[mm]),
							ifelse(sign(pmatrix_oo[mm])==1,(-qnorm(pmatrix_oo[mm])+toadd)^2, 
								qnorm(abs(pmatrix_oo[mm]))+toadd),
							0)
					}
					wZtmp<-pmatrixZ_oo
					wZtmp<-rep(wZtmp,each=(nrepeat*length(treatmentids)))
					wZ<-c(wZ,wZtmp)
				}
			}
			orderwZ<-rep(rep(seq(1:nexp),each=(nrepeat*length(treatmentids))),ntime)
			wZ<-wZ[order(orderwZ,decreasing=FALSE)]
			normtmp[,"TF"]<-normtmp[,ii]
			normtmp[,"wZ"]<-wZ
			nRow=ceiling(-1+sqrt(2+nexp))
			nCol=2+nRow
			titlei<-tf_ids[ii]
			NA0<-apply(data.frame(is.na(normtmp[,"TF"]),normtmp[,"TF"]==0),
				1,function(x){
					rr<-ifelse(any(x==TRUE),TRUE,FALSE)
				return(rr)
			})
			if(all(NA0)){
				next
			}
			summaryData<-matrix(NA,ncol=6,nrow=(ntime*length(treatmentids)))
			idrow=0
			for (aa in 1:ntime){
				for (bb in 1:length(treatmentids)){
					idrow=idrow+1
					tmp<-normtmp[which(normtmp[,"time"]==timevector[aa] & normtmp[,"treatment"]==treatmentids[bb]),]
					if(median(tmp[which(tmp[,"TF"]!=0),"TF"],na.rm=TRUE)=="NaN"){
						mean_w=NA
					} else {
						TF<-tmp[!is.na(tmp[,"TF"]),"TF"]
						wZ<-tmp[!is.na(tmp[,"TF"]),"wZ"]
						zz<-TF*wZ
						mean_w<-ifelse(sum(wZ)!=0,sum(zz)/sum(wZ),NA)
					}
					summaryData[idrow,1]=mean_w
					nrep=ifelse(sum(!is.na(tmp[which(tmp[,"TF"]!=0),"TF"]))==0,1,
						sum(!is.na(tmp[which(tmp[,"TF"]!=0),"TF"])))
					summaryData[idrow,4]=nrep
					if(sd(tmp[which(tmp[,"TF"]!=0),"TF"],na.rm=TRUE)=="NaN"){
						sd_w=NA
					} else {
						TF<-tmp[!is.na(tmp[,"TF"]),"TF"]
						wZ<-tmp[!is.na(tmp[,"TF"]),"wZ"]
						zz<-TF*wZ
						mean_w<-ifelse(sum(wZ)!=0,sum(zz)/sum(wZ),NA)
						var_w<-ifelse(sum(wZ)!=0,(nrep*sum(wZ*((TF-mean_w)^2)))/((nrep-1)*sum(wZ)),NA)
						nonzerowZ<-sum(which(wZ>0))
						var_w<-ifelse(sum(wZ)!=0,
							(sum(wZ*((TF-mean_w)^2))/((nonzerowZ-1)*sum(wZ)/nonzerowZ)),NA)
						sd_w=var_w^0.5
					}
					summaryData[idrow,2]=sd_w
					summaryData[idrow,3]=sd_w/nonzerowZ^0.5
					summaryData[idrow,5]=mean_w-1.96*(sd_w/nonzerowZ^0.5)
					summaryData[idrow,6]=mean_w+1.96*(sd_w/nonzerowZ^0.5)
				}
			}
			summaryData<-data.frame(time=rep(timevector,each=length(treatmentids)),treatment=rep(treatmentids,ntime),summaryData)
			colnames(summaryData)<-c("time","treatment","TF","sd_w","se","nrep","lcl","ucl")
			if (ii==1){
				summaryDataWeighted<-summaryData
			} else {
				summaryDataWeighted<-rbind(summaryDataWeighted,summaryData)
			}
			if (all(!is.finite(summaryData[,"TF"]))){
				next
			}
			pp<-ggplot(data=summaryData,
				aes(x= time,y= TF,colour= treatment,group=treatment))+
				geom_line(na.rm=TRUE)+
				geom_smooth(data=summaryData,aes(ymin = lcl, ymax = ucl,
					fill=treatment), stat="identity")+
				xlab("") +
				ylab("") +
				labs(title =titlei)+
				theme(panel.background = element_rect(fill="transparent"))+
				theme(strip.background=element_rect(fill="white"))+
				theme(axis.ticks=element_line(colour="black"))+
				theme(axis.title=element_text(size=16))+
				theme(axis.text=element_text(size=16,colour="black",face="bold"))+
				theme(axis.title=element_text(size=16))+
				theme(strip.text=element_text(size=16,face="bold"))+
				theme(legend.position="none")+
				theme(panel.border=element_rect(colour="black",fill="transparent"))
			if(is.null(valueColors)){
				pp<-pp+scale_colour_brewer(type="qual",palette=1,name = "Treatment")
			} else {
				pp<-pp+scale_fill_manual(name = "Treatment", values = valueColors)+
					scale_colour_manual(name = "Treatment", values = valueColors)
			}
			assign(paste("plot_",tf_ids[ii],sep=""),pp)
		}
		
		idtf=0
		count=0
		while (idtf<nTF){
			pdf(file = paste(folder.output,"/","Part_",count,"_",pdfname,sep=""),
				paper="a4r",width=11,height=8)
			grid.newpage()
			pushViewport(viewport(layout = grid.layout(4, 5)))
			for (ll in 1:4){
				for (mm in 1:5){
					idtf<-idtf+1
					try(print(get(paste("plot_",tf_ids[idtf],sep="")), 
						vp = vplayout(ll, mm)),TRUE)
				}
			}
			dev.off()
			count=count+1
		}
	} 
	summaryDataWeighted<-cbind(rep(tf_ids,each=(ntime*length(treatmentids))),summaryDataWeighted)
	write.table(summaryDataWeighted, file=paste(folder.output,"/summaryDataWeighted.txt",sep=""), sep="\t",
		row.names = FALSE,col.names = TRUE)	
}

# plotBoxPlot
#============================================================================================

PlotBoxPlot<-function(dataAll=dataFluc,inputparameters=inputParameters,
	expCond=expCond,cex.main=2,cex.lab=2,cex.axis=2){
  folder.output=inputparameters$folder.output
  library (RColorBrewer)
  norm<-dataAll$dataScaleFree
  pdf(file = paste(folder.output,"/Boxplot.pdf",sep=""),paper="a4r",width=11,height=8)
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  norm[,1:nTF]<-apply(norm[,1:nTF],c(1,2),function(x) {r<-ifelse(x==0,NA,x)})
  mycol<-c(brewer.pal(11,"RdYlGn"),brewer.pal(11,"RdYlBu"),brewer.pal(11,"PRGn"))[ntreatment]
  for(i in 1:nTF){
	norm[,"TF"]<-norm[,i]
      if (all(is.na(norm[,"TF"]))){
		next
	}
      ylim<-c(min(norm[,"TF"],na.rm=TRUE),max(norm[,"TF"],na.rm=TRUE))
	nRow=ceiling(-1+sqrt(2+nexp))
   	nCol=2+nRow
      plot.new()
	par(mfrow=c(nRow,nCol),oma = c(0, 0, 3, 0),mar=c(3,7,3,1))
      for (j in 1:nexp){
            dataj<-norm[which(norm[,"experiment"]==j),]
        	boxplot(TF~treatment,data=dataj,
			col=mycol,ylab="Normalized Intensity",main=paste("Experiment ",j,sep=""),
			cex.axis=cex.axis,cex.lab=cex.lab,cex.main=cex.main,ylim=ylim)
	}
	mtext(tf_ids[i], outer = TRUE, cex = 1.5)
   }
   dev.off()
}

# plotPairs
#============================================================================================

PlotPairs<-function(dataAll=dataFluc,inputparameters=inputParameters,
	expCond=expCond,cex.main=2,cex.lab=2,cex.axis=2){
  folder.output=inputparameters$folder.output
  library (RColorBrewer)
  norm<-dataAll$dataScaleFree
  ntreatment<-inputparameters$ntreat
  tf_ids<-expCond$tfNames
  nTF<-inputparameters$nTF
  nexp<-inputparameters$nexp
  nrepeat<-inputparameters$nrepeat
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  norm[,1:nTF]<-apply(norm[,1:nTF],c(1,2),function(x) {r<-ifelse(x==0,NA,x)})
  pdf(file = paste(folder.output,"/TFPairsplot.pdf",sep=""),paper="a4r",width=11,height=8)
  mycol<-brewer.pal(ntreatment,"Blues")
  for(i in 1:nTF){
      par(oma=c(0,3,0,0))
	norm[,"TF"]<-norm[,i]
      if (all(is.na(norm[,"TF"]))){
		next
	}
	panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
	{
    		usr <- par("usr"); on.exit(par(usr))
    		par(usr = c(0, 1, 0, 1))
    		r <- cor(x, y,use="pairwise.complete.obs")
    		txt <- format(c(r, 0.123456789), digits=digits)[1]
    		txt <- paste(prefix, txt, sep="")
    		if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    		text(0.5, 0.5, txt, cex = cex.cor * r)
	}
	panel.line<-function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    	cex = 1, col.line = "red") 
	{
    		points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    		ok <- is.finite(x) & is.finite(y)
    		if (any(ok)) 
        		abline(a=0,b=1,col = col.line)
	}
      tmp1<-norm[,c(tf_ids[i],"experiment")]
	tmp2<-data.frame(unstack(tmp1),norm[which(norm[,"experiment"]==1),"treatment"])
	colnames(tmp2)<-c(paste("Experiment",seq(1,nexp,by=1),sep=" "),"treatment")
      if (nexp>1){
		tmp2[,1:nexp]<-apply(tmp2[,1:nexp],c(1,2),function(x){ifelse(x==0,NA,x)})
      	tmp3<-tmp2[,1:nexp]
      	tmp3<-tmp3[,apply(tmp3,2,function(x){!(all(is.na(x)))})]
		if (!is.null(dim(tmp3))){
			pairs(tmp3,col=mycol[unclass(tmp2[,"treatment"])],
				upper.panel=panel.cor,lower.panel=panel.line,main=tf_ids[i])
      	}
	}    
   }
   dev.off()
   pdf(file = paste(folder.output,"/Pairsplot.pdf",sep=""),paper="a4r",width=11,height=8)
   mycol<-brewer.pal(ntreatment,"Blues")
   library(plyr)
   test<-data.frame(stack(norm[,1:nTF]),
	rep(norm[,"time"],nTF),rep(norm[,"treatment"],nTF),
	rep(norm[,"experiment"]))
   colnames(test)<-c("Act","TF","time","treatment","experiment")
   tmp0<-ddply(test,.(experiment,time,treatment,TF),summarise, 
	mean=mean(Act,na.rm=TRUE))
   for(i in 1:ntime){
      par(oma=c(0,3,0,0))
	panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
	{
    		usr <- par("usr"); on.exit(par(usr))
    		par(usr = c(0, 1, 0, 1))
    		r <- cor(x, y,use="pairwise.complete.obs")
    		txt <- format(c(r, 0.123456789), digits=digits)[1]
    		txt <- paste(prefix, txt, sep="")
    		if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    		text(0.5, 0.5, txt, cex = cex.cor * r)
	}
	panel.line<-function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    	cex = 1, col.line = "red") 
	{
    		points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    		ok <- is.finite(x) & is.finite(y)
    		if (any(ok)) 
        		abline(a=0,b=1,col = col.line)
	}
	tmp2<-data.frame(unstack(tmp0[which(tmp0[,"time"]==timevector[i]),
		c("mean","experiment")]),
		tmp0[which(tmp0[,"time"]==timevector[1]),"treatment"])
	colnames(tmp2)<-c(paste("Exp",seq(1,nexp,by=1),sep=" "),"treatment")
      if (nexp>1){
      	tmp3<-tmp2[,1:nexp]
      	tmp3<-tmp3[,apply(tmp3,2,function(x){!(all(is.na(x)))})]
		if (!is.null(dim(tmp3))){
			pairs(tmp3,col=mycol[unclass(tmp2[,"treatment"])],
				upper.panel=panel.cor,lower.panel=panel.line,
				main=paste("Time",timevector[i],sep=""))
		}
     }
   }
   dev.off()
}

# tableOrganization
#============================================================================================

  tableOrganization<-function(dataFluc=dataAllFluc,
	inputparameters=inputParameters,expCond=expCond){
  #Inputs for the function
	folder.output=inputparameters$folder.output
	tfControl<-inputparameters$tfControl
	tfBlank<-inputparameters$tfBlank
	tfPc<-inputparameters$tfPc
	tfIds<-expCond$tfNames
	treatIds<-expCond$treatmentNames
	normTableInd<-dataFluc$dataScaleFree
	summaryTable<-dataFluc$summaryTable
	summaryBcgTable<-dataFluc$summaryBcgTable
	ntreat<-inputparameters$ntreat
	nTF<-inputparameters$nTF
	normTable<-normTableInd[,1:nTF]
	TF2remove<-inputparameters$TF2remove
	TFExp2removeCond<-inputparameters$TFExp2remove
	nrepeatave<-inputparameters$nrepeatave
	# Generate a text file with the experiments to remove
	TFExpSummary<-ddply(normTableInd,.(experiment,treatment,time),function(x){
		uu<-apply(x[,1:nTF],2,function(y){
				rr<-y[!is.na(y)]
				if (length(rr)==0){
					ss<-NA
				} else {
					tt<-rr[rr!=0]
					ss<-ifelse((length(tt)>(0.5*nrepeatave)),TRUE,FALSE)
				}
				return(ss)
			})
		return(uu)
	})
	TFExp2removefile<-ddply(TFExpSummary,.(experiment),function(x){
		zz<-apply(x[4:ncol(x)],2,function(y){
			if (all(!is.na(y)) && all(y==TRUE)){
				vv<-TRUE
			} else {
				vv<-ifelse(all(is.na(y)),NA,FALSE)
			}
			return(vv)
			})
		return(zz)})
	TFExp2removefilefinal<-NULL	
	for (i in 1:nTF){
		if (tfIds[i]==tfBlank){
			next
		}
		tempexp2rem<-TFExp2removefile[which(TFExp2removefile[,(i+1)]==FALSE),1]
		TFExp2removefilefinal<-rbind(TFExp2removefilefinal,
			cbind(rep(tfIds[i],length(tempexp2rem)),tempexp2rem))
	}
	TFExp2removefilefinal<-data.frame(TFExp2removefilefinal)
	write.table(TFExp2removefilefinal,file="TFExp2remove.txt",sep="\t",
		col.names = FALSE,row.names=FALSE)
	# Remove all the combinations of TF and experiments that should be removed
	summaryTable2<-summaryTable
	summaryBcgTable2<-summaryBcgTable
	if (!is.null(TFExp2removeCond)){
		if (TFExp2removeCond==TRUE){
			TFExp2remove<-TFExp2removefilefinal
		} else {
			TFExp2remove<-read.table(TFExp2removeCond,sep="\t")
		}
		for (jj in 1:nrow(TFExp2remove)){
			normTableInd[which(normTableInd$experiment==TFExp2remove[jj,2]),
				grep(TFExp2remove[jj,1],colnames(normTableInd))]<-NA
			summaryTable2[which(summaryTable2$Experiment==TFExp2remove[jj,2]),
				grep(TFExp2remove[jj,1],colnames(summaryTable2))]<-0
			summaryBcgTable2[which(summaryBcgTable2$Experiment==TFExp2remove[jj,2]),
				grep(TFExp2remove[jj,1],colnames(summaryBcgTable2))]<-0
		}
	}
	TFpresent<-apply(normTableInd[,1:nTF],2,function(x){
		rr<-TRUE
		if (all(is.na(x))){
			rr<-FALSE
		}
		if (all(x[!is.na(x)]==0)){
			rr<-FALSE
		}
		return(rr)
	})
	TFpresent[which(names(TFpresent)==tfControl)]<-FALSE
	TFpresent[which(names(TFpresent)==tfBlank)]<-FALSE
	TFpresent[match(TF2remove,TFpresent)]<-FALSE
	normTable2<-normTableInd[,1:nTF][,which(TFpresent==TRUE)]
	normTableInd2<-data.frame(normTable2,normTableInd[,(nTF+1):ncol(normTableInd)])
	posid<-grep(tfIds[1],colnames(summaryTable2))
	summaryTable3<-summaryTable2[,posid:ncol(summaryTable2)][,which(TFpresent==TRUE)]
	summaryBcgTable3<-summaryTable2[,posid:ncol(summaryBcgTable2)][,which(TFpresent==TRUE)]
	colnames(normTableInd2)<-c(tfIds[TFpresent],colnames(normTableInd)[(nTF+1):ncol(normTableInd)])
	# Identify all the repeats that are not present in any treatment and time
	temprep<-tapply(
		apply(normTable2,1,function(x){
			sum(x,na.rm=TRUE)
		}),normTableInd2[,"repeat"],function(y){
				!all(y==0)
	})
	rowstokeep<-!is.na(match(normTableInd2[,"repeat"],rep(1:nrepeat)[temprep]))
	if (all(rowstokeep==FALSE)){
		print("No data satisfy all the normalization criteria")
		stop
	}
	normTable2<-normTable2[rowstokeep,]
	normTableInd2<-normTableInd2[rowstokeep,]
	# Identify those experiments that should be removed
	exprm<-c()
	for (ee in 1:nexp){
		ini<-grep("TA_above_bcg",colnames(summaryTable))
		fin<-grep("PC_above_bcg",colnames(summaryTable))
		tempexpCT<-apply(summaryTable[which(
			summaryTable[,"Experiment"]==ee),ini:fin],1,function(x){
				above<-as.numeric(x[1])
				r<-ifelse ((above==0 || is.na(above)), TRUE,FALSE)
				return(r)
			}) 
		tempexpPC<-apply(summaryTable[which(
			summaryTable[,"Experiment"]==ee),ini:fin],1,function(x){
				above<-as.numeric(x[2])
				r<-ifelse ((above==0 || is.na(above)), TRUE,FALSE)
				return(r)
			})
		# Avoid the criteria of PC in case there is none
		if (is.na(tfPc)){
			tempexpPC<-FALSE
		}
		# An experiment is removed if any of the controls are not above
		# the background or if all the PC are below the background
		# for all the conditions and times of a given experiment
		if (any(tempexpCT)||all(tempexpPC)){
			tempexp<-TRUE
		} else {
			tempexp<-FALSE
		} 
		exprm<-c(exprm,tempexp)
	}
	# Identify the experiments to remove because all the data are NA or below
	# the background (signified as 0)
	tempexp<-tapply(apply(normTable2,1,function(x){
			sum(x,na.rm=TRUE)}),
			normTableInd2[,"experiment"],function(y){
			rr<-ifelse(all(y==0),TRUE,FALSE)
	})
	exprm<-apply(cbind(exprm,tempexp),1,function(x){
		rr<-ifelse(any(x==TRUE),TRUE,FALSE)
		return(rr)
	})
	if (is.null(exp2remove) && any(exprm==TRUE)){
		exp2remove<-rep(1:nexp)[exprm]
	} 
	if (!is.null(exp2remove) && any(exprm==TRUE)){
		exp2remove<-unique(c(exp2remove,rep(1:nexp)[exprm]))
		exp2remove<-exp2remove[order(exp2remove,decreasing=FALSE)]
	}
	if (length(exp2remove)>0){
		rowstokeep<-is.na(match(normTableInd2[,"experiment"],exp2remove))
		normTable2<-normTable2[rowstokeep,]
		normTableInd2<-normTableInd2[rowstokeep,]	
		# remove the same experiments from the summaryTable and from the 
		# summaryBcgTable
		rowstokeep<-is.na(match(summaryTable2[,"Experiment"],exp2remove))
		summaryTable4<-summaryTable3[rowstokeep,]
		summaryTable2[!rowstokeep,posid:ncol(summaryTable2)]<-0
		rowstokeep<-is.na(match(summaryBcgTable2[,"Experiment"],exp2remove))
		summaryBcgTable4<-summaryBcgTable3[rowstokeep,]
		summaryBcgTable2[!rowstokeep,posid:ncol(summaryBcgTable2)]<-0
	}
	dataFluc[["normTable2"]]<-normTable2
	dataFluc[["normTableInd2"]]<-normTableInd2
  # Update and order the bookkeeping tables. It identifies the
  # TFs that are present for each treatment and experiment at all the
  # times
  orderTable<-matrix(NA,nrow=ntreat,ncol=nTF)
  orderBcgTable<-matrix(NA,nrow=ntreat,ncol=nTF)
  orderBcGArray<-array(NA,dim=c(ntreat,nexp,nTF))
  for (rr in 1:ntreat){
	tempsum<-summaryTable2[which(summaryTable2[,"Treatment"]==treatIds[rr]),]
	tempbcgsum<-summaryBcgTable2[which(summaryBcgTable2[,"Treatment"]==treatIds[rr]),]
	for (ss in 1:nTF){
		temp<-tapply(tempsum[,which(colnames(tempsum)==paste(tfIds[ss],"in_experiment",sep="_"))],
			tempsum[,"Experiment"],function(x){
				if (all(x==0)){
					reSum<-NA
				} else {
					r<-sum(as.numeric(x),na.rm=TRUE)
					reSum<-ifelse(r>0,TRUE,FALSE) 
				}
				return(reSum)
				}
		)
	    tempbcg<-tapply(tempbcgsum[,which(colnames(tempbcgsum)==paste(tfIds[ss],"aboveBcg",sep="_"))],
				tempbcgsum[,"Experiment"],function(x){
					if (all(is.na(x))){
						reSum<-NA
					} else {
						r<-sum(as.numeric(x),na.rm=TRUE)
						reSum<-ifelse(r>0,TRUE,FALSE)
					}
                    return(reSum)
					}
		)
		orderTable[rr,ss]<-sum(temp,na.rm=NA)
		orderBcgTable[rr,ss]<-sum(tempbcg,na.rm=NA)
		orderBcGArray[rr,,ss]<-tempbcg
	}
  }
  # Identify in which specific experiments we can consider each TF to be present (it was in the
  # array and it is above the background and the corresponding experiment has not been excluded
  tfexppresent<-list()
  for (ii in 1: nTF){
    if (nexp==1){
		temptfii<-seq(1,nexp,by=1)
		ordertf<-rep(FALSE,nexp)
		if (all(is.na(orderBcGArray[,1,ii]))){
			ordertf[i]<-FALSE
		} else {
			tt<-sum(orderBcGArray[,1,ii],na.rm=TRUE)
			ordertf[i]<-ifelse(tt>0,TRUE,FALSE)
		}			
	} else {
		temptfii<-seq(1,nexp,by=1)
		ordertf<-apply(orderBcGArray[,,ii],2,function(x){
			if (all(is.na(x))){
				rr<-FALSE
			} else {
				tt<-sum(x,na.rm=TRUE)
				rr<-ifelse(tt>0,TRUE,FALSE)
			}
		})
		temptfii<-temptfii[ordertf]	
    }
	tfexppresent[[tfIds[ii]]]<-temptfii
  }
  colnames(orderTable)<-tfIds
  rownames(orderTable)<-treatIds
  colnames(orderBcgTable)<-tfIds
  rownames(orderBcgTable)<-treatIds
  averageexp<-apply(orderBcgTable,2,min)
  remainTF<-tfIds[which(averageexp>=minexp & TFpresent==TRUE)]
  TFpaf<-apply(normTable2,2,function(x){
	rr<-ifelse(all(is.na(x)),FALSE,TRUE)
	return(rr)
  })
  TFtoexclude<-na.omit(match(colnames(normTable2)[!TFpaf],remainTF))
  if (length(TFtoexclude)>0){
	remainTF<-remainTF[-TFtoexclude]
  }
  dataFluc[["exp2remove"]]<-exp2remove
  dataFluc[["orderTable"]]<-orderTable
  dataFluc[["orderBcgTable"]]<-orderBcgTable
  dataFluc[["orderBcGArray"]]<-orderBcGArray
  ipos<-grep("repeat",colnames(normTableInd2))
  normTableInd3<-normTableInd2[,c(match(remainTF,colnames(normTableInd2)),ipos:ncol(normTableInd2))]
  dataFluc[["normTableInd3"]]<-normTableInd3
  dataFluc[["tfexppresent"]]<-tfexppresent
  dataFluc[["tfIds0"]]<-expCond$tfNames
  dataFluc[["tfIds"]]<-remainTF
  dataFluc[["nrepeat0"]]<-inputparameters$nrepeat
  dataFluc[["nrepeat"]]<-max(normTableInd3[,"repeat"],na.rm=TRUE)
  dataFluc[["nTF0"]]<-inputparameters$nTF
  dataFluc[["nTF"]]<-length(remainTF)
  write.table(normTableInd3, file = paste(folder.output,"/PriorStats.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = NA, qmethod = c("escape", "double"))
  write.table(orderTable, file = paste(folder.output,"/orderTable.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = NA, qmethod = c("escape", "double"))
  write.table(orderBcgTable, file = paste(folder.output,"/orderBcgTable.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = NA, qmethod = c("escape", "double"))
  return(dataFluc)
}

# chol2
#============================================================================================

chol2<-function(x,pivot=FALSE,LINPACK=pivot,tol=-1){
    if (is.complex(x)) 
        stop("complex matrices not permitted at present")
    else if (!is.numeric(x)) 
        stop("non-numeric argument to 'chol'")
    if (is.matrix(x)) {
        if (nrow(x) != ncol(x)) 
            stop("non-square matrix in 'chol'")
        n <- nrow(x)
    }
    else {
        if (length(x) != 1L) 
            stop("non-matrix argument to 'chol'")
        n <- 1L
    }
    n <- as.integer(n)
    if (is.na(n)) 
        stop("invalid nrow(x)")
    if (!is.double(x)) 
        storage.mode(x) <- "double"
    if (pivot) {
        xx <- x
        xx[lower.tri(xx)] <- 0
        z <- .Fortran("dchdc", x = xx, n, n, double(n), piv = integer(n), 
            as.integer(pivot), rank = integer(1L), DUP = FALSE, 
            PACKAGE = "base")
        if (z$rank < n) 
            if (!pivot) 
                stop("matrix not positive definite")
            else warning("matrix not positive definite")
        robj <- z$x
        if (pivot) {
            attr(robj, "pivot") <- z$piv
            attr(robj, "rank") <- z$rank
            if (!is.null(cn <- colnames(x))) 
                colnames(robj) <- cn[z$piv]
        }
        robj
    }
    else {
        z <- .Fortran("chol", x = x, n, n, v = matrix(0, nrow = n, 
            ncol = n), info = integer(1L), DUP = FALSE, PACKAGE = "base")
        if (z$info) 
            stop("non-positive definite matrix in 'chol'")
        z$v
    }
}

# LimmaStat
#============================================================================================

LimmaStat<-function(table6=table6,table7=table7,repeats=repeats,timepoint=timepoint,
	treatment=treatment,array=array,nFactor=nFactor,nameUntreated=nameUntreated,ntime=ntime,
	nrepeat=nrepeatave,ntreat=ntreat,treatmentNames=treatmentNames,timeVector=timeVector,nTF=nTF,
	narray=narray,tfIds=tfIds,treatmentnamesind=treatmentnamesind,timevectorind=timevectorind,
	blockexp=blockexp,method="fdr",customcontrast=customcontrast,treattmp=treattmp){
  sampleIds<-rep(rep(seq(1,ntime*ntreat,by=1),each=nrepeat),narray)
  design<-model.matrix(~0+factor(sampleIds))
  namesGen<-data.frame(Treatment=treatmentnamesind,Time=as.character(timevectorind))
  namesDesign<-apply(namesGen,1,function(x){paste(x[1],"_T",x[2],sep="")})
  colnames(design) <- namesDesign[array==1 & repeats==1]
  dataMatrix<-as.matrix(t(table6))
  #dataMatrix<-apply(dataMatrix,c(1,2),function(x){ifelse(x==0,NA,x)})
  if (!blockexp || narray==1){
	fit<- lmFit(dataMatrix, design=design)
  } else{
	block <- rep(seq(1,narray,by=1),each=ntime*nrepeat*ntreat)
	dupcor <- duplicateCorrelation(dataMatrix,design,block=block)
	fit <- lmFit(dataMatrix,design=design,block=block,correlation=dupcor$consensus)
  }
  fit<- eBayes(fit)
  contrasts<-NULL
  for (jj in 1:(ntreat-1)){
	for (ii in (jj+1):ntreat){
		contrasts = c(contrasts,paste(namesDesign[repeats == 1 & array == 1 & treatment == jj], "-", 
			namesDesign[repeats == 1 & array == 1 & treatment == ii], sep = ""))
	}
  }
  contrastsTime = paste(namesDesign[repeats == 1 & array == 1 & timepoint !=1], "-", 
	namesDesign[repeats == 1 & array == 1 & timepoint != ntime], sep = "")
  contrastsInitTime = paste(namesDesign[repeats == 1 & array == 1 & timepoint !=1], "-", 
	rep(namesDesign[repeats == 1 & array == 1 & timepoint== 1],(ntime-1)),sep = "")
  contrasts<-makeContrasts(contrasts=contrasts, levels=design)
  contrastsTime<-makeContrasts(contrasts=contrastsTime,levels=design)
  contrastsInitTime<-makeContrasts(contrasts=contrastsInitTime,levels=design)
  fit2<-eBayes(contrasts.fit(fit, contrasts))
  fitTime<-eBayes(contrasts.fit(fit, contrastsTime))
  fitInitTime<-eBayes(contrasts.fit(fit, contrastsInitTime))
  pvalue<-apply(fit2$p.value,2,function(x){p.adjust(x,method=method)})
  pvalue<-pvalue*apply(sign(fit2$coefficient),c(1,2),function(x){ifelse(x==0,1,x)})
  pvalueTime<-apply(fitTime$p.value,2,function(x){p.adjust(x,method=method)})
  pvalueTime<-pvalueTime*apply(sign(fitTime$coefficient),c(1,2),function(x){ifelse(x==0,1,x)})
  pvalueInitTime<-apply(fitInitTime$p.value,2,function(x){p.adjust(x,method=method)})
  pvalueInitTime<-pvalueInitTime*apply(sign(fitInitTime$coefficient),c(1,2),
	function(x){ifelse(x==0,1,x)})
  if (any(nameUntreated!="none")){ 
	treattoinclude<-unique(treatmentnamesind)[-match(unique(treattmp),unique(treatmentnamesind))]
	treattoinclude<-treatmentNames[match(treattoinclude,treatmentNames)]
	sampleControlIds<- rep(rep(seq(1,(ntime*length(treattoinclude)),by=1),each=nrepeat),narray)
	designControl <- model.matrix(~ 0+factor(sampleControlIds))
	namesGenvsControl<-data.frame(Treatment=treattmp,Time=as.character(timevectorind))
	namesDesignvsControl<-apply(namesGenvsControl,1,function(x){paste(x[1],"_T",x[2],sep="")})
	namesGenControl<-data.frame(Treatment=treatmentnamesind[which(is.na(match(treatmentnamesind,treattmp)))],
		Time=as.character(timevectorind[which(is.na(match(treatmentnamesind,treattmp)))]))
    namesDesignControl<-apply(namesGenControl,1,function(x){paste(x[1],"_T",x[2],sep="")})
	arrayControl<-array[which(is.na(match(treatmentnamesind,treattmp)))]
	repeatsControl<-repeats[which(is.na(match(treatmentnamesind,treattmp)))]
	timepointControl<-timepoint[which(is.na(match(treatmentnamesind,treattmp)))]
	colnames(designControl) <- namesDesignControl[arrayControl==1 & repeatsControl==1]
	dataMatrixvsControl<-as.matrix(t(table7))[,which(is.na(match(treatmentnamesind,treattmp)))]
	dataMatrixControl<-as.matrix(t(table7))[,which(is.na(match(treatmentnamesind,treattmp)))]
	if (!blockexp || narray==1){
		fitvsControl<-lmFit(dataMatrix, design=design)
		fitControl <- lmFit(dataMatrixControl, design=designControl)
	} else{
		block <- rep(seq(1,narray,by=1),each=ntime*nrepeat*ntreat)
		dupcor <- duplicateCorrelation(dataMatrix,design,block=block)
		fitvsControl <- lmFit(dataMatrix,design=design,block=block,correlation=dupcor$consensus)
		blockControl <- rep(seq(1,narray,by=1),each=ntime*nrepeat*length(treattoinclude))
		dupcorControl <- duplicateCorrelation(dataMatrixControl,designControl,block=blockControl)
		fitControl <- lmFit(dataMatrixControl,design=designControl,block=blockControl,
			correlation=dupcor$consensus)
	}
	contrastsvsControl<-paste(namesDesignControl[repeatsControl == 1 & arrayControl == 1],
		"-", namesDesignvsControl[which(treatmentnamesind!=treattmp)][repeatsControl == 1 & arrayControl == 1], sep = "")
	contrastsControlTime = paste(namesDesignControl[repeatsControl == 1 & arrayControl == 1 & timepointControl !=1], "-", 
		namesDesignControl[repeatsControl == 1 & arrayControl == 1 & timepointControl != ntime], sep = "")
	contrastsControlInitTime = paste(namesDesignControl[repeatsControl == 1 & arrayControl == 1 & timepointControl !=1], "-", 
		rep(namesDesignControl[repeatsControl == 1 & arrayControl == 1 & timepointControl== 1],(ntime-1)),sep = "")
	contrastsvsControl<-makeContrasts(contrasts=contrastsvsControl,levels=design)
	contrastsControlTime<-makeContrasts(contrasts=contrastsControlTime,levels=designControl)
	contrastsControlInitTime<-makeContrasts(contrasts=contrastsControlInitTime,levels=designControl)
	fitvsControl2<-eBayes(contrasts.fit(fitvsControl, contrastsvsControl))
	fitControlTime<-eBayes(contrasts.fit(fitControl, contrastsControlTime))
	fitControlInitTime<-eBayes(contrasts.fit(fitControl, contrastsControlInitTime))
	pvaluevsControl<-apply(fitvsControl2$p.value,2,function(x){p.adjust(x,method=method)})
    pvaluevsControl<-pvaluevsControl*apply(sign(fitvsControl2$coefficient),
		c(1,2),function(x){ifelse(x==0,1,x)})
	pvalueControlTime<-apply(fitControlTime$p.value,2,function(x){p.adjust(x,method=method)})
    pvalueControlTime<-pvalueControlTime*apply(sign(fitControlTime$coefficient),
		c(1,2),function(x){ifelse(x==0,1,x)})
	pvalueControlInitTime<-apply(fitControlInitTime$p.value,2,function(x){p.adjust(x,method=method)})
    pvalueControlInitTime<-pvalueControlInitTime*apply(sign(fitControlInitTime$coefficient),
		c(1,2),function(x){ifelse(x==0,1,x)})
  } else {
    pvaluevsControl<-c()
	pvalueControlTime<-c()
	pvalueControlInitTime<-c()
  }
  statsLimma<-list(pvalue=pvalue,pvalueTime=pvalueTime,pvalueInitTime=pvalueInitTime,
	pvaluevsControl=pvaluevsControl,pvalueControlTime=pvalueControlTime,
	pvalueControlInitTime=pvalueControlInitTime)
  return(statsLimma)
}
# WrappingLimmaStat
#============================================================================================
  WrapLimmaStat<-function(table6=table6,table7=table7,repeats=repeats,timepoint=timepoint,treatment=treatment,
		array=array,nFactor=nFactor,nameUntreated=nameUntreated,ntime=ntime,nrepeat=nrepeatave,ntreat=ntreat,
		treatmentNames=treatmentNames,timeVector=timeVector,nTF=nTF,narray=narray,tfIds=tfIds,
		treatmentnamesind=treatmentnamesind,timevectorind=timevectorind,
		blockexp=blockexp,method=method,customcontrast=customcontrast,treattmp=treattmp,
		folder.output=folder.output){
  # Call Limma Stat to run obtain p-values for each contrast
  LimmaResults<-LimmaStat(table6=table6,table7=table7,repeats=repeats,timepoint=timepoint,treatment=treatment,
	array=array,nFactor=nFactor,nameUntreated=nameUntreated,ntime=ntime,nrepeat=nrepeatave,ntreat=ntreat,treatmentNames=treatmentNames,
	timeVector=timeVector,nTF=nTF,narray=narray,tfIds=tfIds,treatmentnamesind=treatmentnamesind,timevectorind=timevectorind,
	blockexp=blockexp,method=method,customcontrast=customcontrast,treattmp=treattmp)
  # Write results
  write.table(LimmaResults[["pvalue"]], file =paste(folder.output,"/TreatmentComparisons.txt", sep=""),
	append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names =TRUE, 
	col.names = TRUE, qmethod = c("escape", "double"))
  pvalueNLCA<-apply(LimmaResults[["pvalue"]],c(1,2),function(x){
	if (abs(x)>pvaluelimma||is.na(x)){			
		r<-0
	} else {
		r<-1*sign(x)
	}
	return(r)
  })
  write.table(pvalueNLCA, file = paste(folder.output,"/TreatmentComparisonsNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
  write.table(LimmaResults[["pvalueTime"]], file = paste(folder.output,"/TimeComparisons.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
  pvalueTimeNLCA<-apply(LimmaResults[["pvalueTime"]],c(1,2),function(x){
	if (abs(x)>pvaluelimma||is.na(x)){
		r<-0
	} else {
		r<-1*sign(x)
	}
	return(r)
  })
  pvalueTimeNLCA<-data.frame(pvalue=rep(1,nrow(pvalueTimeNLCA)),pvalueTimeNLCA)
  write.table(pvalueTimeNLCA, file = paste(folder.output,"/TimeComparisonsNLCA.txt", sep=""),append = FALSE, 
	quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names =TRUE, 
	col.names = TRUE, qmethod = c("escape", "double"))
  write.table(LimmaResults[["pvalueInitTime"]], file = paste(folder.output,"/InitialTimeComparison.txt",sep=""), append = FALSE, 
	quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, 
	qmethod = c("escape", "double"))
  pvalueInitTimeNLCA<-apply(LimmaResults[["pvalueInitTime"]],c(1,2),function(x){
	if (abs(x)>pvaluelimma||is.na(x)){
		r<-0
	} else {
		r<-1*sign(x)
	}
	return(r)
  })
  pvalueInitTimeNLCA<-data.frame(pvalue=rep(1,nrow(pvalueInitTimeNLCA)),pvalueInitTimeNLCA)
  write.table(pvalueInitTimeNLCA, file = paste(folder.output,"/InitialTimeComparisonsNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
# In case controls are present in the experiments
  if(any(nameUntreated!="none")){
	write.table(LimmaResults[["pvaluevsControl"]], file = paste(folder.output,"/ComparisonsVsControl.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
	pvaluevsControlNLCA<-apply(LimmaResults[["pvaluevsControl"]],c(1,2),function(x){
		if (abs(x)>pvaluelimma||is.na(x)){
			r<-0
		} else {
			r<-1*sign(x)
		}
		return(r)
	})
	pvaluevsControlNLCA<-data.frame(pvalue=rep(1,nrow(pvaluevsControlNLCA)),pvaluevsControlNLCA)
	write.table(pvaluevsControlNLCA, file = paste(folder.output,"/ComparisonsVsControlNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
	write.table(LimmaResults[["pvalueControlTime"]], file = paste(folder.output,"/TimeComparisonsVsControl.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
	pvalueControlTimeNLCA<-apply(LimmaResults[["pvalueControlTime"]],c(1,2),function(x){
		if (abs(x)>pvaluelimma||is.na(x)){
			r<-0
		} else {
			r<-1*sign(x)
		}
		return(r)
	})
	pvalueControlTimeNLCA<-data.frame(pvalue=rep(1,nrow(pvalueControlTimeNLCA)),pvalueControlTimeNLCA)
	write.table(pvalueControlTimeNLCA, file = paste(folder.output,"/TimeComparisonsVsControlNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))
	write.table(LimmaResults[["pvalueControlInitTime"]], file =paste(folder.output,"/InitialTimeComparisonVsControl.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
	pvalueControlInitTimeNLCA<-apply(LimmaResults[["pvalueControlInitTime"]],c(1,2),function(x){
		if (abs(x)>pvaluelimma||is.na(x)){
			r<-0
		} else {
			r<-1*sign(x)
	    }
		return(r)
	})
	pvalueControlInitTimeNLCA<-data.frame(pvalue=rep(1,nrow(pvalueControlInitTimeNLCA)),pvalueControlInitTimeNLCA)
	write.table(pvalueControlInitTimeNLCA, file = paste(folder.output,"/InitialTimeComparisonsVsControlNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
			eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))
  } else {
	pvaluevsControlNLCA<-NULL
	pvalueControlTimeNLCA<-NULL
	pvalueControlInitTimeNLCA<-NULL
  }
  LimmaResults[["pvalueNLCA"]]=pvalueNLCA
  LimmaResults[["pvalueTimeNLCA"]]=pvalueTimeNLCA
  LimmaResults[["pvalueInitTimeNLCA"]]=pvalueInitTimeNLCA
  LimmaResults[["pvaluevsControlNLCA"]]=pvaluevsControlNLCA
  LimmaResults[["pvalueControlTimeNLCA"]]=pvalueControlTimeNLCA
  LimmaResults[["pvalueControlInitTimeNLCA"]]=pvalueControlInitTimeNLCA	  
  return(LimmaResults)
 }

# WrappingLimma
#============================================================================================
WrappingLimma<-function(dataFluc=dataFluc,inputparameters=inputparameters,
	expCond=expCond,method="fdr"){
  folder.output=inputparameters$folder.output
  library('limma')
  require(limma)
  normTable<-dataFluc$normTableInd3
  tfIds<-dataFluc$tfIds
  nTF<-dataFluc$nTF
  minexp<-inputparameters$minexp
  narray<-minexp
  ntreat<-inputparameters$ntreat
  nrepeat0<-dataFluc$nrepeat
  nrepeatave<-inputparameters$nrepeatave
  ntime<-inputparameters$ntime
  nrowmin<-nrepeat0*ntime*ntreat*minexp
  nameUntreated<-inputparameters$nameUntreated
  orderBcgTable<-dataFluc$orderBcgTable[,match(
	dataFluc$tfIds,colnames(dataFluc$orderBcgTable))]
  nexpprs<-apply(orderBcgTable,2,min)
  meanexp<-mean(nexpprs,na.rm=TRUE)
  treatmentNames<-expCond$treatmentNames
  timeVector<-inputparameters$timeVector
  nFactor<-inputparameters$nFactor
  blockexp<-inputparameters$blockexp
  methodlimma<-inputparameters$methodlimma
  pvaluelimma<-inputparameters$pvaluelimma
  customcontrast<-inputparameters$customcontrast
  metanalysis<-inputparameters$metanalysis
  pvaluemetanalysis<-inputparameters$pvaluemetanalysis
  # Generate initial tables
  explist<-unique(normTable$experiment)
  nexp<-length(explist)
  table6<-dataFluc$dataNLCA[,1:nTF]
  repeats<-rep(seq(1,nrepeatave,by=1),ntime*ntreat*nexp)
  treatment<-rep(rep(seq(1,ntreat,by=1),each=nrepeatave),ntime*nexp)
  treatmentnamesind<-rep(rep(treatmentNames,each=nrepeatave),ntime*nexp)
  timepoint<-rep(rep(seq(1,ntime,by=1),each=nrepeatave*ntreat),nexp)
  timevectorind<-rep(rep(timeVector,each=nrepeatave*ntreat),nexp)
  array<-rep(explist,each=ntime*ntreat*nrepeatave)
  totallen<-ntime*ntreat*nrepeatave*nexp
  factorIds<-matrix(NA,nrow=nFactor,ncol=totallen)
  for (ii in 1:nFactor){
	nfactorii<-length(expCond$levelsNames[[ii]])
	factorIds[ii,]<-rep(seq(1,nfactorii,by=1),each=(totallen/nfactorii))
  }
  namesGen<-data.frame(Experiment=array,Treatment=treatmentnamesind,Time=as.character(timevectorind),Repeat=repeats)
  rownames(table6)<-apply(namesGen,1,function(x){paste("Exp",x[1],"_",x[2],"_T",x[3],"_",x[4],sep="")})
  table6updated<-table6
  table6updated<-apply(table6,c(1,2),function(x){
	rr<-ifelse(is.finite(x),x,NA)
	return(rr)
  })
  table6updated[which(table6updated=="NaN")]<-NA
  table6ind<-data.frame(table6updated,repeats=repeats,treatment=treatmentnamesind,time=timevectorind,experiment=array,
	t(factorIds))
  colnames(table6ind)<-c(colnames(normTable)[1:nTF],"repeats","treatment","time","experiment",
	colnames(normTable)[(ncol(normTable)+1-nFactor):ncol(normTable)])
  rownames(table6ind)<-c()
  dataFluc[["table6ind"]]<-table6ind
  if (any(nameUntreated!="none")){
	table7<-matrix(NA,nrow(table6),ncol(table6))
	table6NA<-apply(table6,c(1,2),function(x){ifelse(x==0,NA,x)})
	fid<-which(nameUntreated!="none")
	treattemp<-NULL
	if (nFactor>1){
		for (rr in 1:length(treatmentnamesind)){
			treattemp<-rbind(treattemp,
				unlist(strsplit(treatmentnamesind[rr],"\\.")))
		}
		factortemp<-data.frame(treattemp)
		for (tt in 1:length(fid)){
			treattemp[,fid[tt]]<-nameUntreated[fid]
		}
		treattmp<-apply(treattemp,1,function(x){
			paste(x,collapse=".")}
		)
	} else {
		treattemp<-treatmentnamesind
		factortemp<-treatmentnamesind
		treattemp<-rep(nameUntreated,length(treatmentnamesind))
		treattmp<-treattemp
	}
	for (aa in 1:ntime){
		for (bb in 1:ntreat){
			for (cc in 1:narray){
				for (dd in 1:nTF){
					treattempjj<-unique(treattmp[which(
							as.logical(array == explist[cc] & treatment == bb & timepoint == aa))]
					)
					meantmp<-mean(table6NA[which(as.logical(array == explist[cc] & treatment==which(treatmentNames==treattempjj) & timepoint == aa)),dd],
						na.rm=TRUE)
					meantmp<-ifelse(meantmp=="NaN",0,meantmp)
					table7[which(as.logical(array == explist[cc] & treatment == bb & timepoint == aa)),dd]<-
						table6[which(as.logical(array == explist[cc] & treatment == bb & timepoint == aa)),dd]-meantmp
				}
			}	
		}
	}
	colnames(table7)<-colnames(table6)
	rownames(table7)<-rownames(table6)
	table7updated<-table7
	table7updated[which(table7updated=="NaN")]<-NA
	table7ind<-data.frame(table7updated,repeats,timevectorind,treatmentnamesind,array,
		factortemp)
	colnames(table7ind)<-c(colnames(table7),"repeat","time","treatment","experiment",namesFactors)
	rownames(table7ind)<-c()
	dataFluc[["table7ind"]]<-table7ind
 } else {
	table7<-NULL
	treattmp<-NULL
 }
# Determine if bootstrapping is required
  error<-abs(nexpprs)-meanexp
  if (all(error<.Machine$double.eps ^ 0.5)){
	BootReq<-FALSE
  } else {
	BootReq<-TRUE
  }
  nboot=inputparameters$nboot
  ntreatment<-inputparameters$ntreat
  expCond$tfNames<-dataAllFluc$tfIds
  ntime<-inputparameters$ntime
  timevector<-inputparameters$timeVector
  treatmentids<-expCond$treatmentNames
  inputparameters$nTF<-dataAllFluc$nTF
  inputparameters$nexp<-nexp
  inputparameters$nrepeat<-nrepeatave
  # All the experiments have the same number of TFs and normalization can be applied
  if (!BootReq && !metanalysis){
	LimmaResults<-WrapLimmaStat(table6=table6,table7=table7,repeats=repeats,timepoint=timepoint,treatment=treatment,
		array=array,nFactor=nFactor,nameUntreated=nameUntreated,ntime=ntime,nrepeat=nrepeatave,ntreat=ntreat,
		treatmentNames=treatmentNames,timeVector=timeVector,nTF=nTF,narray=narray,tfIds=tfIds,
		treatmentnamesind=treatmentnamesind,timevectorind=timevectorind,
		blockexp=blockexp,method=method,customcontrast=customcontrast,treattmp=treattmp,
		folder.output=folder.output)
	if(any(nameUntreated!="none")){
		dataFluc[["LimmaResults"]]<-list(pvalue=LimmaResults[["pvalue"]],pvalueTime=LimmaResults[["pvalueTime"]],
			pvalueInitTime=LimmaResults[["pvalueInitTime"]],pvaluevsControl=LimmaResults[["pvaluevsControl"]],
			pvalueControlTime=LimmaResults[["pvalueControlTime"]],pvalueControlInitTime=LimmaResults[["pvalueControlInitTime"]])
	} else {
		dataFluc[["LimmaResults"]]<-list(pvalue=LimmaResults[["pvalue"]],pvalueTime=LimmaResults[["pvalueTime"]],
			pvalueInitTime=LimmaResults[["pvalueInitTime"]])
	}
  }
  # When meta-analysis is required because regular normalization do not work
  if (metanalysis){
	ResMetaAnalysis<-list()
	table6complete<-table6
	table7complete<-table7
	arraycomplete<-array
	nTFcomplete<-nTF
	tfIdscomplete<-tfIds
	narraycomplete<-narray
	repeatscomplete=repeats
	timepointcomplete=timepoint
	treatmentcomplete=treatment
	treatmentnamesindcomplete=treatmentnamesind
	timevectorindcomplete=timevectorind
	blockexpcomplete<-blockexp
	treattmpcomplete<-treattmp
	for (ll in 1:nexp){
		print(ll)
		table6<-table6complete[arraycomplete==explist[ll],]
		NaNPres<-which(colSums(table6)=="NaN")
		NaPres<-which(colSums(table6)=="NA")
		remove_NANaN<-unique(c(NaNPres,NaPres))
		remove_NANaN<-remove_NANaN[order(remove_NANaN)]
		if (length(remove_NANaN)>0){
			table6<-table6[,-remove_NANaN]
		}
		if(is.null(dim(table6))|| all (is.na(table6))){
			next
		}
		if(!is.null(table7)){
			table7<-table7complete[arraycomplete==explist[ll],]
			NaNPres<-which(colSums(table7)=="NaN")
			NaPres<-which(colSums(table7)=="NA")
			remove_NANaN<-unique(c(NaNPres,NaPres))
			remove_NANaN<-remove_NANaN[order(remove_NANaN)]
			if (length(remove_NANaN)>0){
			table7<-table7[,-remove_NANaN]
			}
			if(is.null(dim(table7))|| all (is.na(table6))){
				table7<-NULL
			}
		}
		repeats=repeatscomplete[arraycomplete==explist[ll]]
		timepoint=timepointcomplete[arraycomplete==explist[ll]]
		treatment=treatmentcomplete[arraycomplete==explist[ll]]
		array<-arraycomplete[arraycomplete==explist[ll]]
		array<-rep(1,length(array))
		tfIds<-names(which(colSums(table6)!="NaN"))
		nTF<-length(tfIds)
		treatmentnamesind=treatmentnamesindcomplete[arraycomplete==explist[ll]]
		timevectorind=timevectorindcomplete[arraycomplete==explist[ll]]
		blockexp<-FALSE
		narray=1
		treattmp<-treattmpcomplete[arraycomplete==explist[ll]]
		# Statistics for the first experiment
		LimmaResults<-WrapLimmaStat(table6=table6,table7=table7,repeats=repeats,timepoint=timepoint,treatment=treatment,
			array=array,nFactor=nFactor,nameUntreated=nameUntreated,ntime=ntime,nrepeat=nrepeatave,ntreat=ntreat,
			treatmentNames=treatmentNames,timeVector=timeVector,nTF=nTF,narray=narray,tfIds=tfIds,
			treatmentnamesind=treatmentnamesind,timevectorind=timevectorind,
			blockexp=blockexp,method=method,customcontrast=customcontrast,treattmp=treattmp,
			folder.output=folder.output)
		tmppvalue<-as.data.frame(LimmaResults[["pvalue"]])
		tmppvalueTime<-as.data.frame(LimmaResults[["pvalueTime"]])
		tmppvalueInitTime<-as.data.frame(LimmaResults[["pvalueInitTime"]])	
		tmppvaluevsControl<-as.data.frame(LimmaResults[["pvaluevsControl"]])
		tmppvalueControlTime<-as.data.frame(LimmaResults[["pvalueControlTime"]])		
		tmppvalueControlInitTime<-as.data.frame(LimmaResults[["pvalueControlInitTime"]])				
		if(ll==1){
				pvalueMA<-array(NA,dim=c(nTFcomplete,ncol(tmppvalue),nexp))
				pvalueTimeMA<-array(NA,dim=c(nTFcomplete,ncol(tmppvalueTime),nexp))
				pvalueInitTimeMA<-array(NA,dim=c(nTFcomplete,ncol(tmppvalueInitTime),nexp))
				pvaluevsControlMA<-array(NA,dim=c(nTFcomplete,ncol(tmppvaluevsControl),nexp))
				pvalueControlTimeMA<-array(NA,dim=c(nTFcomplete,ncol(tmppvalueControlTime),nexp))
				pvalueControlInitTimeMA<-array(NA,dim=c(nTFcomplete,ncol(tmppvalueControlInitTime),nexp))
		}
		pvalueMA[match(rownames(tmppvalue),tfIdscomplete),,ll]<-as.matrix(tmppvalue)
		pvalueTimeMA[match(rownames(tmppvalueTime),tfIdscomplete),,ll]<-as.matrix(tmppvalueTime)
		pvalueInitTimeMA[match(rownames(tmppvalueInitTime),tfIdscomplete),,ll]<-as.matrix(tmppvalueInitTime)
		pvaluevsControlMA[match(rownames(tmppvaluevsControl),tfIdscomplete),,ll]<-as.matrix(tmppvaluevsControl)
		pvalueControlTimeMA[match(rownames(tmppvalueControlTime),tfIdscomplete),,ll]<-
			as.matrix(tmppvalueControlTime)
		pvalueControlInitTimeMA[match(rownames(tmppvalueControlInitTime),tfIdscomplete),,ll]<-
			as.matrix(tmppvalueControlInitTime)		
	}
	#Calculate the chi-square test.
	# The function is modified to account for the sign of the differences that are being compared
	FischerSq<-function(p=p,pMAcutoff=pvaluemetanalysis){
		if (all(is.na(p))){
			p.val<-NA
		} else {
			if (length(which(p>=0))>0){
				ppos<-p[which(p>=0)]
			} else {
				ppos<-1
			}
			if (length(which(p<0))>0){
				pneg<-p[which(p<0)]
			} else {
				pneg<-1
			}
			total_df<-length(ppos)+length(pneg)
			Xsqpos<-(-2*sum(log(ppos)))
			p.valpos<-1-pchisq(Xsqpos, df = 2*total_df)
			Xsqneg<-(-2*sum(log(abs(pneg))))
			p.valneg<-1-pchisq(Xsqneg, df = 2*total_df)
			if(p.valpos<=pMAcutoff && p.valneg<=pMAcutoff){
				p.val<-1
			}	
			if(p.valpos<=pMAcutoff && p.valneg>pMAcutoff){
				p.val<-p.valpos
			}	
			if(p.valpos>pMAcutoff && p.valneg<=pMAcutoff){
				p.val<-(-p.valneg)
			}
			if(p.valpos>pMAcutoff && p.valneg>pMAcutoff){
				p.val<-ifelse(p.valpos>p.valneg,(-p.valneg),p.valpos)
			}
		}
		return(p.val)
	}
	zTransform<-function(p=p,pMAcutoff=pvaluemetanalysis){
		if (all(is.na(p))){
			p.val<-NA
		} else {
			pfinite<-p[is.finite(p)]
			npfinite<-length(pfinite)
			ztrans<-data.frame(
				p.value=pfinite,
				signp=sign(pfinite)
			)
			ztrans["Zvalue"]=apply(ztrans,1,function(x){
				rr<-qnorm(abs(x[1]))
				return(rr)
			})
			Zscorepos=sum(ztrans[which(ztrans$p.value>0),"Zvalue"])/(npfinite^0.5)
			Zscoreneg=sum(ztrans[which(ztrans$p.value<0),"Zvalue"])/(npfinite^0.5)
			p.valpos=pnorm(Zscorepos)
			p.valneg=pnorm(Zscoreneg)
			if(p.valpos<=pMAcutoff && p.valneg<=pMAcutoff){
				p.val<-1
			}	
			if(p.valpos<=pMAcutoff && p.valneg>pMAcutoff){
				p.val<-p.valpos
			}	
			if(p.valpos>pMAcutoff && p.valneg<=pMAcutoff){
				p.val<-(-p.valneg)
			}
			if(p.valpos>pMAcutoff && p.valneg>pMAcutoff){
				p.val<-ifelse(p.valpos>p.valneg,(-p.valneg),p.valpos)
			}
		}
		return(p.val)
	}
	pvalueMA_entirearray<-aperm(pvalueMA,c(2,1,3))
	pvalueMA_entirearray<-aperm(pvalueMA_entirearray,c(1,3,2))
	dim(pvalueMA_entirearray)<-c((dim(pvalueMA)[2]*dim(pvalueMA)[3]),dim(pvalueMA)[1])
	colnames(pvalueMA_entirearray)<-tfIdscomplete
	expids<-unique(table6ind[,"experiment"])
	pvalueMA_entirearray<-data.frame(experiment=rep(expids,each=dim(pvalueMA)[2]),comparisons=rep(colnames(tmppvalue),length(expids)),pvalueMA_entirearray)
	write.table(pvalueMA_entirearray, file = paste(folder.output,"/TreatmentComparisonsMA_entirearray.txt",sep=""), 
		append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", 
		dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
	pvalueFinal<-apply(pvalueMA,c(1,2),function(x){zTransform(p=x)})
	colnames(pvalueFinal)<-colnames(tmppvalue)
	rownames(pvalueFinal)<-tfIdscomplete
	write.table(pvalueFinal, file = paste(folder.output,"/TreatmentComparisonsMA.txt",sep=""), 
		append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", 
		dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
	pvalueTimeFinal<-apply(pvalueTimeMA,c(1,2),function(x){zTransform(p=x)})
	colnames(pvalueTimeFinal)<-colnames(tmppvalueTime)
	rownames(pvalueTimeFinal)<-tfIdscomplete
	write.table(pvalueTimeFinal, file = paste(folder.output,"/TimeComparisonsMA.txt", sep=""),
		append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, 
		qmethod = c("escape", "double"))
	pvalueInitTimeMA_entirearray<-aperm(pvalueInitTimeMA,c(2,1,3))
	pvalueInitTimeMA_entirearray<-aperm(pvalueInitTimeMA_entirearray,c(1,3,2))
	dim(pvalueInitTimeMA_entirearray)<-c((dim(pvalueInitTimeMA)[2]*dim(pvalueInitTimeMA)[3]),dim(pvalueInitTimeMA)[1])
	colnames(pvalueInitTimeMA_entirearray)<-tfIdscomplete
	pvalueInitTimeMA_entirearray<-data.frame(experiment=rep(expids,each=dim(pvalueInitTimeMA)[2]),comparisons=rep(colnames(tmppvalueInitTime),length(expids)),pvalueInitTimeMA_entirearray)
	write.table(pvalueInitTimeMA_entirearray, file = paste(folder.output,"/InitialTimeComparisonsMA_entirearray.txt",sep=""), 
		append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", 
		dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
	pvalueInitTimeFinal<-apply(pvalueInitTimeMA,c(1,2),function(x){zTransform(p=x)})
	colnames(pvalueInitTimeFinal)<-colnames(tmppvalueInitTime)
	rownames(pvalueInitTimeFinal)<-tfIdscomplete
	write.table(pvalueInitTimeFinal, file = paste(folder.output,"/InitialTimeComparisonsMA.txt",sep=""),
		append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
		row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
	if (any(nameUntreated!="none")){
		pvaluevsControlFinal<-apply(pvaluevsControlMA,c(1,2),function(x){zTransform(p=x)})
		colnames(pvaluevsControlFinal)<-colnames(tmppvaluevsControl)
		rownames(pvaluevsControlFinal)<-tfIdscomplete
		write.table(pvaluevsControlFinal, file = paste(folder.output,"/ComparisonsVsControl.txt", sep=""),
			append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
			row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
		pvalueControlTimeFinal<-apply(pvalueControlTimeMA,c(1,2),function(x){zTransform(p=x)})
		colnames(pvalueControlTimeFinal)<-colnames(tmppvalueControlTime)
		rownames(pvalueControlTimeFinal)<-tfIdscomplete
		write.table(pvalueControlTimeFinal, file = paste(folder.output,"/TimeComparisonsVsControl.txt", sep=""),
			append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
			row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
		pvalueControlInitTimeFinal<-apply(pvalueControlInitTimeMA,c(1,2),function(x){zTransform(p=x)})
		colnames(pvalueControlInitTimeFinal)<-colnames(tmppvalueControlInitTime)
		rownames(pvalueControlInitTimeFinal)<-tfIdscomplete
		write.table(pvalueControlTimeFinal, file = paste(folder.output,"/InitTimeComparisonsVsControl.txt", sep=""),
			append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
			row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
	} else {
		pvaluevsControlFinal<-NULL
		pvalueControlTimeFinal<-NULL
		pvalueControlInitTimeFinal<-NULL
	}
	dataFluc[["LimmaResults"]]<-list(pvalue=pvalueFinal,pvalueTime=pvalueTimeFinal,
		pvalueInitTime=pvalueInitTimeFinal,pvaluevsControl=pvaluevsControlFinal,
		pvalueControlTime=pvalueControlTimeFinal,pvalueControlInitTime=pvalueControlInitTimeFinal,
		pvalueMA=pvalueMA,pvalueTimeMA=pvalueTimeMA,
		pvalueInitTimeMA=pvalueInitTimeMA,pvaluevsControlMA=pvaluevsControlMA,
		pvalueControlTimeMA=pvalueControlTimeMA,pvalueControlInitTimeMA=pvalueControlInitTimeMA)
  }
  # When bootstrapping is required and there is no need for a meta-analysis study  
  if (BootReq && !metanalysis){
	for (kk in 1:nboot){
		print (kk)
		table6temp<-matrix(NA,nrow=nrowmin,ncol=nTF)
		for (ii in 1:nTF){
			nexpii<-length(dataFluc$tfexppre[[tfIds[ii]]])
			if (nexpii==1){
				expii<-dataFluc$tfexppre[[tfIds[ii]]]
			} else {
				expii<-sample(dataFluc$tfexppre[[tfIds[ii]]],minexp,replace=FALSE)
			}
			table6temp[,ii]<-normTable[!is.na(match(normTable$experiment, expii)),ii]
		}
		tablestack<-as.numeric(table6temp)
		table8<-data.frame(tablestack=tablestack,block=rep(seq(1,length(tablestack)/nrepeat0,by=1),each=nrepeat0))
		table9<-tapply(table8$tablestack,table8$block,function(x){
			gooddata<-x[which(is.finite(x)>0)]
			ngooddata<-length(gooddata)
			if (ngooddata==nrepeatave){
				r<-gooddata
			} else {
				if (ngooddata==0){
					r<-x[1:nrepeatave]
				}
				if (ngooddata<nrepeatave){ #missing values are substitute by the mean
					meannona<-mean(gooddata)
					r<-c(gooddata, rep(meannona,(nrepeatave-ngooddata)))
				} 
				if (sum(!is.na(x))>nrepeatave){ # randomly remove the excess of data
					r<-sample(gooddata,nrepeatave,replace=FALSE)
				}
			}
		return(r)
		})
		table6<-matrix(unlist(table9),ncol=nTF)
		colnames(table6)<-tfIds
		repeats<-rep(seq(1,nrepeatave,by=1),ntime*ntreat*narray)
		treatment<-rep(rep(seq(1,ntreat,by=1),each=nrepeatave),ntime*narray)
		treatmentnamesind<-rep(rep(treatmentNames,each=nrepeatave),ntime*narray)
		timepoint<-rep(rep(seq(1,ntime,by=1),each=nrepeatave*ntreat),narray)
		timevectorind<-rep(rep(timeVector,each=nrepeatave*ntreat),narray)
		array<-rep(seq(1,narray,by=1),each=ntime*ntreat*nrepeatave)
		totallen<-ntime*ntreat*nrepeatave*narray
		factors<-matrix(NA,nrow=nFactor,ncol=totallen)
		for (ii in 1:nFactor){
			nfactorii<-length(expCond$levelsNames[[ii]])
			factors[ii,]<-rep(seq(1,nfactorii,by=1),each=(totallen/nfactorii))
		}
		namesGen<-data.frame(Experiment=array,Treatment=treatmentnamesind,Time=as.character(timevectorind))
		rownames(table6)<-apply(namesGen,1,function(x){paste("Exp",x[1],"_",x[2],"_T",x[3],sep="")})
		if (any(nameUntreated!="none")){
			table7<-matrix(NA,nrow(table6),ncol(table6))
			table6NA<-apply(table6,c(1,2),function(x){ifelse(x==0,NA,x)})
			fid<-which(nameUntreated!="none")
			treattemp<-NULL
			if (nFactor>1){
				for (rr in 1:length(treatmentnamesind)){
					treattemp<-rbind(treattemp,
					unlist(strsplit(treatmentnamesind[rr],"\\.")))
				}
				factortemp<-data.frame(treattemp)
				for (tt in 1:length(fid)){
					treattemp[,fid[tt]]<-nameUntreated[fid]
				}
				treattmp<-apply(treattemp,1,function(x){
					paste(x,collapse=".")}
				)
			} else {
				treattemp<-treatmentnamesind
				factortemp<-treatmentnamesind
				treattemp<-rep(nameUntreated,length(treatmentnamesind))
				treattmp<-treattemp
			}
			for (aa in 1:ntime){
				for (bb in 1:ntreat){
					for (cc in 1:narray){
						for (dd in 1:nTF){
							treattempjj<-unique(treattmp[which(
								as.logical(array == cc & treatment == bb & timepoint == aa))]
							)
							meantmp<-mean(table6NA[which(as.logical(array == cc & treatment==which(treatmentNames==treattempjj) & timepoint == aa)),dd],
								na.rm=TRUE)
							meantmp<-ifelse(meantmp=="NaN",0,meantmp)
							table7[which(as.logical(array == cc & treatment == bb & timepoint == aa)),dd]<-
								table6[which(as.logical(array == cc & treatment == bb & timepoint == aa)),dd]-meantmp
						}
					}	
				}
			}
			colnames(table7)<-colnames(table6)
			rownames(table7)<-rownames(table6)
			table7ind<-data.frame(table7,repeats,timevectorind,treatmentnamesind,array,
				factortemp)
			colnames(table7ind)<-c(colnames(table7),"repeat","time","treatment","experiment",namesFactors)
		} else {
			treattmp<-NULL
		}
		LimmaResults<-LimmaStat(table6=table6,table7=table7,repeats=repeats,timepoint=timepoint,treatment=treatment,
			array=array,nFactor=nFactor,nameUntreated=nameUntreated,ntime=ntime,nrepeat=nrepeatave,ntreat=ntreat,treatmentNames=treatmentNames,
			timeVector=timeVector,nTF=nTF,narray=narray,tfIds=tfIds,treatmentnamesind=treatmentnamesind,timevectorind=timevectorind,
			blockexp=blockexp,method=methodlimma,customcontrast=customcontrast,treattmp=treattmp)
		if (kk==1){
			LimmaBootstrap<-LimmaResults
			LimmaBootstrap[["pvalue"]]<-as.numeric(LimmaResults[["pvalue"]])
			LimmaBootstrap[["pvalueTime"]]<-as.numeric(LimmaResults[["pvalueTime"]])
			LimmaBootstrap[["pvalueInitTime"]]<-as.numeric(LimmaResults[["pvalueInitTime"]])
			if (is.null(LimmaResults[["pvaluevsControl"]])){
				LimmaBootstrap[["pvaluevsControl"]]<-c()
			} else {
				LimmaBootstrap[["pvaluevsControl"]]<-as.numeric(LimmaResults[["pvaluevsControl"]])
			}
			if (is.null(LimmaResults[["pvalueControlTime"]])){
				LimmaBootstrap[["pvalueControlTime"]]<-c()
			} else {
				LimmaBootstrap[["pvalueControlTime"]]<-as.numeric(LimmaResults[["pvalueControlTime"]])
			}
			if (is.null(LimmaResults[["pvalueControlInitTime"]])){	
				LimmaBootstrap[["pvalueControlInitTime"]]<-c()
			} else {
				LimmaBootstrap[["pvalueControlInitTime"]]<-as.numeric(LimmaResults[["pvalueControlInitTime"]])
			}
		} else {
			LimmaBootstrap[["pvalue"]]<-cbind(LimmaBootstrap[["pvalue"]],
				as.numeric(LimmaResults[["pvalue"]]))
			LimmaBootstrap[["pvalueTime"]]<-cbind(LimmaBootstrap[["pvalueTime"]],
				as.numeric(LimmaResults[["pvalueTime"]]))
			LimmaBootstrap[["pvalueInitTime"]]<-cbind(LimmaBootstrap[["pvalueInitTime"]],
				as.numeric(LimmaResults[["pvalueInitTime"]]))
			if (is.null(LimmaResults[["pvaluevsControl"]])){
				LimmaBootstrap[["pvaluevsControl"]]<-c()
			} else {
				LimmaBootstrap[["pvaluevsControl"]]<-cbind(LimmaBootstrap[["pvaluevsControl"]],
					as.numeric(LimmaResults[["pvaluevsControl"]]))
			}
			if (is.null(LimmaResults[["pvalueControlTime"]])){
				LimmaBootstrap[["pvalueControlTime"]]<-c()
			} else {
				LimmaBootstrap[["pvalueControlTime"]]<-cbind(LimmaBootstrap[["pvalueControlTime"]],
					as.numeric(LimmaResults[["pvalueControlTime"]]))
			}
			if (is.null(LimmaResults[["pvalueControlInitTime"]])){	
				LimmaBootstrap[["pvalueControlInitTime"]]<-c()
			} else {
				LimmaBootstrap[["pvalueControlInitTime"]]<-cbind(LimmaBootstrap[["pvalueControlInitTime"]],
					as.numeric(LimmaResults[["pvalueControlInitTime"]]))
			}
		}
	}
	LimmaBootstrap[["pvalue"]]<-data.frame(TFs=rep(tfIds,ncol(LimmaResults[["pvalue"]])),
		Contrasts=rep(colnames(LimmaResults[["pvalue"]]),each=nTF),LimmaBootstrap[["pvalue"]],
		average=apply(LimmaBootstrap[["pvalue"]],1,function(x){mean(x,na.rm=TRUE)}),
		Q1=apply(LimmaBootstrap[["pvalue"]],1,function(x){quantile(x,0.01,na.rm=TRUE)}),
		Q5=apply(LimmaBootstrap[["pvalue"]],1,function(x){quantile(x,0.05,na.rm=TRUE)}),
		Q10=apply(LimmaBootstrap[["pvalue"]],1,function(x){quantile(x,0.10,na.rm=TRUE)}),
		Q25=apply(LimmaBootstrap[["pvalue"]],1,function(x){quantile(x,0.25,na.rm=TRUE)}),
		Q50=apply(LimmaBootstrap[["pvalue"]],1,function(x){quantile(x,0.50,na.rm=TRUE)}),
		Q75=apply(LimmaBootstrap[["pvalue"]],1,function(x){quantile(x,0.75,na.rm=TRUE)}),
		Q90=apply(LimmaBootstrap[["pvalue"]],1,function(x){quantile(x,0.90,na.rm=TRUE)}),
		Q95=apply(LimmaBootstrap[["pvalue"]],1,function(x){quantile(x,0.95,na.rm=TRUE)}),
		Q99=apply(LimmaBootstrap[["pvalue"]],1,function(x){quantile(x,0.99,na.rm=TRUE)}),
		SIG=apply(LimmaBootstrap[["pvalue"]],1,function(x){
			r<-sign(max(x,na.rm=TRUE)*min(x,na.rm=TRUE))
			t<-ifelse(r==1,TRUE,FALSE)
			return(t)
		})
	)
	LimmaBootstrap[["pvalueTime"]]<-data.frame(TFs=rep(tfIds,ncol(LimmaResults[["pvalueTime"]])),
		Contrasts=rep(colnames(LimmaResults[["pvalueTime"]]),each=nTF),LimmaBootstrap[["pvalueTime"]],
		average=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){mean(x,na.rm=TRUE)}),
		Q1=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){quantile(x,0.01,na.rm=TRUE)}),
		Q5=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){quantile(x,0.05,na.rm=TRUE)}),
		Q10=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){quantile(x,0.10,na.rm=TRUE)}),
		Q25=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){quantile(x,0.25,na.rm=TRUE)}),
		Q50=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){quantile(x,0.50,na.rm=TRUE)}),
		Q75=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){quantile(x,0.75,na.rm=TRUE)}),
		Q90=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){quantile(x,0.90,na.rm=TRUE)}),
		Q95=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){quantile(x,0.95,na.rm=TRUE)}),
		Q99=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){quantile(x,0.99,na.rm=TRUE)}),
		SIG=apply(LimmaBootstrap[["pvalueTime"]],1,function(x){
			r<-sign(max(x,na.rm=TRUE)*min(x,na.rm=TRUE))
			t<-ifelse(r==1,TRUE,FALSE)
			return(t)
		})
	)
	LimmaBootstrap[["pvalueInitTime"]]<-data.frame(TFs=rep(tfIds,ncol(LimmaResults[["pvalueInitTime"]])),
		Contrasts=rep(colnames(LimmaResults[["pvalueInitTime"]]),each=nTF),LimmaBootstrap[["pvalueInitTime"]],
		average=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){mean(x,na.rm=TRUE)}),
		Q1=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){quantile(x,0.01,na.rm=TRUE)}),
		Q5=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){quantile(x,0.05,na.rm=TRUE)}),
		Q10=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){quantile(x,0.10,na.rm=TRUE)}),
		Q25=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){quantile(x,0.25,na.rm=TRUE)}),
		Q50=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){quantile(x,0.50,na.rm=TRUE)}),
		Q75=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){quantile(x,0.75,na.rm=TRUE)}),
		Q90=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){quantile(x,0.90,na.rm=TRUE)}),
		Q95=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){quantile(x,0.95,na.rm=TRUE)}),
		Q99=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){quantile(x,0.99,na.rm=TRUE)}),
		SIG=apply(LimmaBootstrap[["pvalueInitTime"]],1,function(x){
			r<-sign(max(x,na.rm=TRUE)*min(x,na.rm=TRUE))
			t<-ifelse(r==1,TRUE,FALSE)
			return(t)
		})
	)
	if (is.null(LimmaResults[["pvaluevsControl"]])){
		LimmaBootstrap[["pvaluevsControl"]]<-c()
	} else {
		LimmaBootstrap[["pvaluevsControl"]]<-data.frame(TFs=rep(tfIds,ncol(LimmaResults[["pvaluevsControl"]])),
			Contrasts=rep(colnames(LimmaResults[["pvaluevsControl"]]),each=nTF),LimmaBootstrap[["pvaluevsControl"]],
			average=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){mean(x,na.rm=TRUE)}),
			Q1=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){quantile(x,0.01,na.rm=TRUE)}),
			Q5=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){quantile(x,0.05,na.rm=TRUE)}),
			Q10=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){quantile(x,0.10,na.rm=TRUE)}),
			Q25=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){quantile(x,0.25,na.rm=TRUE)}),
			Q50=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){quantile(x,0.50,na.rm=TRUE)}),
			Q75=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){quantile(x,0.75,na.rm=TRUE)}),
			Q90=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){quantile(x,0.90,na.rm=TRUE)}),
			Q95=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){quantile(x,0.95,na.rm=TRUE)}),
			Q99=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){quantile(x,0.99,na.rm=TRUE)}),
			SIG=apply(LimmaBootstrap[["pvaluevsControl"]],1,function(x){
				r<-sign(max(x,na.rm=TRUE)*min(x,na.rm=TRUE))
				t<-ifelse(r==1,TRUE,FALSE)
				return(t)
			})
		)
	}
	if (is.null(LimmaResults[["pvalueControlTime"]])){
		LimmaBootstrap[["pvalueControlTime"]]<-c()
	} else {
		LimmaBootstrap[["pvalueControlTime"]]<-data.frame(TFs=rep(tfIds,ncol(LimmaResults[["pvalueControlTime"]])),
			Contrasts=rep(colnames(LimmaResults[["pvalueControlTime"]]),each=nTF),LimmaBootstrap[["pvalueControlTime"]],
			average=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){mean(x,na.rm=TRUE)}),
			Q1=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){quantile(x,0.01,na.rm=TRUE)}),
			Q5=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){quantile(x,0.05,na.rm=TRUE)}),
			Q10=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){quantile(x,0.10,na.rm=TRUE)}),
			Q25=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){quantile(x,0.25,na.rm=TRUE)}),
			Q50=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){quantile(x,0.50,na.rm=TRUE)}),
			Q75=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){quantile(x,0.75,na.rm=TRUE)}),
			Q90=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){quantile(x,0.90,na.rm=TRUE)}),
			Q95=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){quantile(x,0.95,na.rm=TRUE)}),
			Q99=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){quantile(x,0.99,na.rm=TRUE)}),
			SIG=apply(LimmaBootstrap[["pvalueControlTime"]],1,function(x){
				r<-sign(max(x,na.rm=TRUE)*min(x,na.rm=TRUE))
				t<-ifelse(r==1,TRUE,FALSE)
				return(t)
			})
		)
	}
	if (is.null(LimmaResults[["pvalueControlInitTime"]])){
		LimmaBootstrap[["pvalueControlInitTime"]]<-c()
	} else {
		LimmaBootstrap[["pvalueControlInitTime"]]<-data.frame(TFs=rep(tfIds,ncol(LimmaResults[["pvalueControlInitTime"]])),
			Contrasts=rep(colnames(LimmaResults[["pvalueControlInitTime"]]),each=nTF),LimmaBootstrap[["pvalueControlInitTime"]],
			average=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){mean(x,na.rm=TRUE)}),
			Q1=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){quantile(x,0.01,na.rm=TRUE)}),
			Q5=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){quantile(x,0.05,na.rm=TRUE)}),
			Q10=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){quantile(x,0.10,na.rm=TRUE)}),
			Q25=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){quantile(x,0.25,na.rm=TRUE)}),
			Q50=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){quantile(x,0.50,na.rm=TRUE)}),
			Q75=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){quantile(x,0.75,na.rm=TRUE)}),
			Q90=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){quantile(x,0.90,na.rm=TRUE)}),
			Q95=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){quantile(x,0.95,na.rm=TRUE)}),
			Q99=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){quantile(x,0.99,na.rm=TRUE)}),
			SIG=apply(LimmaBootstrap[["pvalueControlInitTime"]],1,function(x){
				r<-sign(max(x,na.rm=TRUE)*min(x,na.rm=TRUE))
				t<-ifelse(r==1,TRUE,FALSE)
				return(t)
			})
		)
	}
	write.table(LimmaBootstrap[["pvalue"]][,c("TFs","Contrasts","average","Q1","Q5","Q10","Q25","Q50","Q75","Q90","Q95","Q99","SIG")],
		file = paste(folder.output,"/TreatmentComparisonsAllData.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
	pvalueFinal<-apply(data.frame(LimmaBootstrap[["pvalue"]][,"Q1"],
		LimmaBootstrap[["pvalue"]][,"Q99"],LimmaBootstrap[["pvalue"]][,"SIG"]),1,function(x){
			if (abs(x[1])>abs(x[2])){
					mm<-x[1]
			} else { 
				mm<-x[2]
			}
			if(!x[3]) mm<-1
			return(mm)
		})
	pvalueFinal<-matrix(pvalueFinal,nrow=nTF)
	colnames(pvalueFinal)<-colnames(LimmaResults[["pvalue"]])
	rownames(pvalueFinal)<-rownames(LimmaResults[["pvalue"]])
	write.table(pvalueFinal, file = paste(folder.output,"/TreatmentComparisons.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
	pvalueFinalNLCA<-apply(pvalueFinal,c(1,2),function(x){
		if ((abs(x)>pvaluelimma)||is.na(x)){
			r<-0
		} else {
			r<-1*sign(x)
		}
		return(r)
	})
	pvalueFinalNLCA<-pvalueFinalNLCA[rowSums(abs(pvalueFinalNLCA))>0,]
	pvalueFinalNLCA<-data.frame(pvalue=rep(1,nrow(pvalueFinalNLCA)),pvalueFinalNLCA)
	write.table(pvalueFinalNLCA, file = paste(folder.output,"/TreatmentComparisonsNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))	
	pvalueTimeFinal<-apply(data.frame(LimmaBootstrap[["pvalueTime"]][,"Q1"],
		LimmaBootstrap[["pvalueTime"]][,"Q99"],LimmaBootstrap[["pvalueTime"]][,"SIG"]),1,function(x){
			if (abs(x[1])>abs(x[2])){
					mm<-x[1]
			} else { 
				mm<-x[2]
			}
			if(!x[3]) mm<-1
			return(mm)
		})
	pvalueTimeFinal<-matrix(pvalueTimeFinal,nrow=nTF)
	colnames(pvalueTimeFinal)<-colnames(LimmaResults[["pvalueTime"]])
	rownames(pvalueTimeFinal)<-rownames(LimmaResults[["pvalueTime"]])
	write.table(pvalueTimeFinal, file = paste(folder.output,"/TimeComparisons.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
	write.table(LimmaBootstrap[["pvalueTime"]][,c("TFs","Contrasts","average","Q1","Q5","Q10","Q25","Q50","Q75","Q90","Q95","Q99","SIG")],
		file = paste(folder.output,"/TimeComparisonsAllData.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
	pvalueTimeFinalNLCA<-apply(pvalueTimeFinal,c(1,2),function(x){
		if ((abs(x)>pvaluelimma)||is.na(x)){
			r<-0
		} else {
			r<-1*sign(x)
		}
		return(r)
	})
	pvalueTimeFinalNLCA<-pvalueTimeFinalNLCA[rowSums(abs(pvalueTimeFinalNLCA))>0,]
	if (!is.null(nrow(pvalueTimeFinalNLCA))){
		pvalueTimeFinalNLCA<-data.frame(pvalue=rep(1,nrow(pvalueTimeFinalNLCA)),pvalueTimeFinalNLCA)
		write.table(pvalueTimeFinalNLCA, file =paste(folder.output,"/TimeComparisonsNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
			eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))
	}
	pvalueInitTimeFinal<-apply(data.frame(LimmaBootstrap[["pvalueInitTime"]][,"Q1"],
		LimmaBootstrap[["pvalueInitTime"]][,"Q99"],LimmaBootstrap[["pvalueInitTime"]][,"SIG"]),1,function(x){
			if (abs(x[1])>abs(x[2])){
					mm<-x[1]
			} else { 
				mm<-x[2]
			}
			if(!x[3]) mm<-1
			return(mm)
		})
	pvalueInitTimeFinal<-matrix(pvalueInitTimeFinal,nrow=nTF)
	colnames(pvalueInitTimeFinal)<-colnames(LimmaResults[["pvalueInitTime"]])
	rownames(pvalueInitTimeFinal)<-rownames(LimmaResults[["pvalueInitTime"]])
	write.table(pvalueInitTimeFinal, file = paste(folder.output,"/InitialTimeComparison.txt",sep=""),append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
	write.table(LimmaBootstrap[["pvalueInitTime"]][,c("TFs","Contrasts","average","Q1","Q5","Q10","Q25","Q50","Q75","Q90","Q95","Q99","SIG")],
		file = paste(folder.output,"/InitialTimeComparisonAllData.txt",sep=""),append = FALSE, quote = FALSE, sep = "\t", 
		eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
	pvalueInitTimeFinalNLCA<-apply(pvalueInitTimeFinal,c(1,2),function(x){
		if ((abs(x)>pvaluelimma)||is.na(x)){
			r<-0
		} else {
			r<-1*sign(x)
		}
		return(r)
	})
	pvalueInitTimeFinalNLCA<-pvalueInitTimeFinalNLCA[rowSums(abs(pvalueInitTimeFinalNLCA))>0,]
	if (!is.null(nrow(pvalueInitTimeFinalNLCA))){
		pvalueInitTimeFinalNLCA<-data.frame(pvalue=rep(1,nrow(pvalueInitTimeFinalNLCA)),pvalueInitTimeFinalNLCA)
		write.table(pvalueInitTimeFinalNLCA, file = paste(folder.output,"/InitialTimeComparisonsNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
			eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
		}
	if (any(nameUntreated!="none")){
		pvaluevsControlFinal<-apply(data.frame(LimmaBootstrap[["pvaluevsControl"]][,"Q1"],
			LimmaBootstrap[["pvaluevsControl"]][,"Q99"],LimmaBootstrap[["pvaluevsControl"]][,"SIG"]),1,function(x){
			if (abs(x[1])>abs(x[2])){
					mm<-x[1]
			} else { 
				mm<-x[2]
			}
			if(!x[3]) mm<-1
			return(mm)
		})	
		pvaluevsControlFinal<-matrix(pvaluevsControlFinal,nrow=nTF)
		colnames(pvaluevsControlFinal)<-colnames(LimmaResults[["pvaluevsControl"]])
		rownames(pvaluevsControlFinal)<-rownames(LimmaResults[["pvaluevsControl"]])
		write.table(pvaluevsControlFinal, file = paste(folder.output,"/ComparisonsVsControl.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
			eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
		write.table(LimmaBootstrap[["pvaluevsControl"]][,c("TFs","Contrasts","average","Q1","Q5","Q10","Q25","Q50","Q75","Q90","Q95","Q99","SIG")],
			file = paste(folder.output,"/ComparisonsVsControl.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
			eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
		pvaluevsControlFinalNLCA<-apply(pvaluevsControlFinal,c(1,2),function(x){
			if ((abs(x)>pvaluelimma)||is.na(x)){
				r<-0
			} else {
				r<-1*sign(x)
			}
			return(r)
		})
		pvaluevsControlFinalNLCA<-pvaluevsControlFinalNLCA[rowSums(abs(pvaluevsControlFinalNLCA))>0,]
		if (!is.null(nrow(pvaluevsControlFinalNLCA))){
			pvaluevsControlFinalNLCA<-data.frame(pvalue=rep(1,nrow(pvaluevsControlFinalNLCA)),pvaluevsControlFinalNLCA)
			write.table(pvaluevsControlFinalNLCA, file = paste(folder.output,"/ComparisonsVsControlNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
				eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
		}
		pvalueControlTimeFinal<-apply(data.frame(LimmaBootstrap[["pvalueControlTime"]][,"Q1"],
			LimmaBootstrap[["pvalueControlTime"]][,"Q99"],LimmaBootstrap[["pvalueControlTime"]][,"SIG"]),1,function(x){
			if (abs(x[1])>abs(x[2])){
					mm<-x[1]
			} else { 
				mm<-x[2]
			}
			if(!x[3]) mm<-1
			return(mm)
		})	
		pvalueControlTimeFinal<-matrix(pvalueControlTimeFinal,nrow=nTF)
		colnames(pvalueControlTimeFinal)<-colnames(LimmaResults[["pvalueControlTime"]])
		rownames(pvalueControlTimeFinal)<-rownames(LimmaResults[["pvalueControlTime"]])
		write.table(pvalueControlTimeFinal, file = paste(folder.output,"/TimeComparisonsVsControl.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t", 
			eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
		write.table(LimmaBootstrap[["pvalueControlTime"]][,c("TFs","Contrasts","average","Q1","Q5","Q10","Q25","Q50","Q75","Q90","Q95","Q99","SIG")], 
			file = paste(folder.output,"/TimeComparisonsVsControlAllData.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t", 
			eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
		pvalueControlTimeFinalNLCA<-apply(pvalueControlTimeFinal,c(1,2),function(x){
			if ((abs(x)>pvaluelimma)||is.na(x)){
				r<-0
			} else {
				r<-1*sign(x)
			}
			return(r)
		})
		pvalueControlTimeFinalNLCA<-pvalueControlTimeFinalNLCA[rowSums(abs(pvalueControlTimeFinalNLCA))>0,]
		if (!is.null(nrow(pvalueControlTimeFinalNLCA))){
			pvalueControlTimeFinalNLCA<-data.frame(pvalue=rep(1,nrow(pvalueControlTimeFinalNLCA)),pvalueControlTimeFinalNLCA)
			write.table(pvalueControlTimeFinalNLCA, file = paste(folder.output,"/TimeComparisonsVsControlNLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
				eol = "\n", na = "NA", dec = ".", row.names =TRUE, col.names = TRUE, qmethod = c("escape", "double"))
		}
		pvalueControlInitTimeFinal<-apply(data.frame(LimmaBootstrap[["pvalueControlInitTime"]][,"Q1"],
			LimmaBootstrap[["pvalueControlInitTime"]][,"Q99"],LimmaBootstrap[["pvalueControlInitTime"]][,"SIG"]),1,function(x){
				if (abs(x[1])>abs(x[2])){
					mm<-x[1]
				} else { 
					mm<-x[2]
				}
				if(!x[3]) mm<-1
				return(mm)
			})	
		pvalueControlInitTimeFinal<-matrix(pvalueControlInitTimeFinal,nrow=nTF)
		colnames(pvalueControlInitTimeFinal)<-colnames(LimmaResults[["pvalueControlInitTime"]])
		rownames(pvalueControlInitTimeFinal)<-rownames(LimmaResults[["pvalueControlInitTime"]])
		write.table(pvalueControlInitTimeFinal, file = paste(folder.output,"/InitialTimeComparisonVsControl.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
			eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
		write.table(LimmaBootstrap[["pvalueControlInitTime"]][,c("TFs","Contrasts","average","Q1","Q5","Q10","Q25","Q50","Q75","Q90","Q95","Q99","SIG")], 
			file = paste(folder.output,"/InitialTimeComparisonVsControlAllData.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
			eol = "\n", na = "NA", dec = ".", row.names =TRUE , col.names = TRUE, qmethod = c("escape", "double"))
		pvalueControlInitTimeFinalNLCA<-apply(pvalueControlInitTimeFinal,c(1,2),function(x){
			if ((abs(x)>pvaluelimma)||is.na(x)){
				r<-0
			} else {
				r<-1*sign(x)
			}
			return(r)
			})
		pvalueControlTimeFinalNLCA<-pvalueControlTimeFinalNLCA[rowSums(abs(pvalueControlTimeFinalNLCA))>0,]
			if (!is.null(nrow(pvalueControlTimeFinalNLCA))){
				pvalueControlTimeFinalNLCA<-data.frame(pvalue=rep(1,nrow(pvalueControlInitTimeFinalNLCA)),pvalueControlInitTimeFinalNLCA)
				write.table(pvalueControlInitTimeFinalNLCA, file = paste(folder.output,"/InitialTimeComparisonVsControlNLCA.txt",sep=""),
					append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names =TRUE, 
					col.names = TRUE, qmethod = c("escape", "double"))
			}
	} else {
			pvaluevsControlFinal<-c()
			pvalueControlTimeFinal<-c()
			pvalueControlInitTimeFinal<-c()
	}
	dataFluc[["LimmaResults"]]<-list(pvalue=pvalueFinal,pvalueTime=pvalueTimeFinal,
		pvalueInitTime=pvalueInitTimeFinal,pvaluevsControl=pvaluevsControlFinal,
		pvalueControlTime=pvalueControlTimeFinal,pvalueControlInitTime=pvalueControlInitTimeFinal,
		LimmaBootstrap=LimmaBootstrap)
  }
  # Plotting only the final experiments that are going to be used in limma. This should be done a posteriori, after limma
  # so you get the p-value and translate them into Z-values for the plotting. That averages are more representative of what it is going on
  plotnorm(dataAll=table6ind,inputparameters=inputparameters,
	expCond=expCond,pdfname="FinalTF.pdf",ytitle="Normalized Intensity")
  if (metanalysis && nFactor==1){
		pmatrix<-pvalueMA
		pmatrixZ<-pvalueFinal
		namespmatrix<-colnames(pvalueFinal)
		plotnormWeigthedAverage(dataAll=table6ind ,inputparameters=inputparameters,
			expCond=expCond,pdfname="FinalTFAverage.pdf",ytitle="Average Normalized Intensity",
			pmatrix=pmatrix,pmatrixZ=pmatrixZ,namespmatrix=namespmatrix)
		plotnormWeigthedAverageSmallArrays(dataAll=table6ind ,inputparameters=inputparameters,
			expCond=expCond,pdfname="FinalTFAverageAllTFtogether.pdf",ytitle="Average Normalized Intensity",
			pmatrix=pmatrix,pmatrixZ=pmatrixZ,namespmatrix=namespmatrix)
  } else{
	plotnormAverage(dataAll=table6ind ,inputparameters=inputparameters,
		expCond=expCond,pdfname="FinalTFAverage.pdf",ytitle="Average Normalized Intensity")
  }
  if (nTF<=20){
	if (!metanalysis && nFactor==1){
		plotnormAverageSmallArrays(dataAll=table6ind ,inputparameters=inputparameters,
			expCond=expCond,pdfname="FinalTFAverageAllTFtogether.pdf",
			ytitle="Average Normalized Intensity")
	}
  }
  if (any(inputparameters$nameUntreated!="none")){
	plotnorm(dataAll=table7ind,
		inputparameters=inputparameters,expCond=expCond,
		pdfname="FinalTFvsControl.pdf",
		ytitle="Normalized Intensity versus Control")
	if (metanalysis && nFactor==1){
		pmatrix<-pvaluevsControlMA
		pmatrixZ<-pvaluevsControlFinal
		namespmatrix<-colnames(pvaluevsControlFinal)
		plotnormWeigthedAverage(dataAll=table7ind ,inputparameters=inputparameters,
			expCond=expCond,pdfname="FinalTFAveragevsControl.pdf",
			ytitle="Average Normalized Intensity versus Control",pmatrix=pmatrix,
			pmatrixZ=pmatrixZ,namespmatrix=namespmatrix)
		plotnormWeigthedAverageSmallArrays(dataAll=table6ind ,inputparameters=inputparameters,
			expCond=expCond,pdfname="FinalTFAverageAllTFtogether.pdf",ytitle="Average Normalized Intensity",
			pmatrix=pmatrix,pmatrixZ=pmatrixZ,namespmatrix=namespmatrix)
	} else{
		plotnormAverage(dataAll=table7ind,
			inputparameters=inputparameters,expCond=expCond,
			pdfname="FinalTFAveragevsControl.pdf",
			ytitle="Average Normalized Intensity versus Control")
	}
	if (nTF<=20){
		if (!metanalysis && nFactor==1){
			plotnormAverageSmallArrays(dataAll=table7ind,
				inputparameters=inputparameters,expCond=expCond,
				pdfname="FinalTFAverageAllTFtogethervsControl.pdf",
				ytitle="Average Normalized Intensity versus Control")
		}
	}
  }
  return(dataFluc)
}

# statAnalysis
#============================================================================================

statAnalysis<-function(dataFluc=dataAllFluc,
	inputparameters=inputParameters,expCond=expCond){
  folder.output=inputparameters$folder.output
  method=inputparameters$methodlimma
  pvalue=inputparameters$pvaluelimma
  treatmentIds=expCond$treatmentNames
  nrepeatave=inputparameters$nrepeatave
  nrepeats=length(unique(dataFluc$normTableInd3[,"repeat"]))
  ntime=length(unique(dataFluc$normTableInd3[,"time"]))
  nexp=length(unique(dataFluc$normTableInd3$experiment))
  nTF=dataFluc$nTF
  ntreat=length(treatmentIds)
  nFactor=inputparameters$nFactor
  namesFactors=inputparameters$namesFactors
  levelNames=expCond$levelNames
  block<-rep(seq(1:(ntime*ntreat*nexp)),each=nrepeats)
  dataTableRaw<-data.frame(dataFluc$normTableInd3, block=block)
  # Table for NLCA. Replace the missing data with the average of the data.
  # If data is in excess, it randomly picked some values
  dataTable<-apply(dataTableRaw[,1:nTF],2,function(x){
		tapply(x,dataTableRaw$block,function(y){
			rr<-y
			gooddata<-y[which(is.finite(y)>0)]
			ngooddata<-length(gooddata)
			if (ngooddata==nrepeatave){
				rr[1:nrepeatave]<-gooddata
			} else {
				if (ngooddata<nrepeatave){ #missing values are substitute by the mean
					meannona<-mean(gooddata)
					rr[1:nrepeatave]<-c(gooddata, rep(meannona,(nrepeatave-ngooddata)))
				} 
				if (sum(!is.na(y))>nrepeatave){ # randomly remove the excess of data
					rr[1:nrepeatave]<-sample(gooddata,nrepeatave,replace=FALSE)
				}
			}
		return(rr)
        })
  })
  dataTable<-data.frame(matrix(unlist(dataTable),ncol=nTF,byrow=0),dataTableRaw[,(nTF+1):(ncol(dataTableRaw)-1)])
  colnames(dataTable)<-c(dataFluc$tfIds,"repeats","time","plate","treatment","experiment",expCond$factorNames)
  dataTable<-dataTable[which(dataTable$repeats<=nrepeatave),]
  dataTable<-data.frame(dataTable[1:nTF],dataTable$repeats,dataTable$time,dataTable$treatment,dataTable[,expCond$factorNames],dataTable$experiment)
  colnames(dataTable)<-c(dataFluc$tfIds,"repeats","time","treatment",expCond$factorNames,"experiment")
  dataTable$time<-rep(rep(seq(1,ntime,by=1),each=(nrepeatave*ntreat)),nexp)
  dataTable$treatment<-rep(rep(seq(1,ntreat,by=1),each=nrepeatave),nexp*ntime)
  write.table(dataTable, file = paste(folder.output,"/DATATABLENLCA.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", 
	eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))
  dataFluc[["dataNLCA"]]<-dataTable
  dataFluc<-WrappingLimma(dataFluc=dataFluc,
  	inputparameters=inputparameters,expCond=expCond,
	method=method)
}




