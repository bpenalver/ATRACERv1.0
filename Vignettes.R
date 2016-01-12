
################################################################################################################
################################################################################################################

############                          ANALYZING LIVING CELL ARRAYS                      ########################
################################################################################################################
################################################################################################################


# This code will take the raw data from a living cell array (LCA) experiment and rearranges them, normalizes them
# and it will provide a list of the most relevant TFs for your study

# This program is a generalization and improvement of the
# original program created and developed in R by Beatriz Penalver Bernabe ##

# Original Code Published In:

# Weiss MS, Penalver Bernabe B, Bellis AD, Broadbelt LJ, Jeruss JS, Shea LD.
# PLoS One. 2010 Nov 17;5(11):e14026

# Bellis AD, Penalver Bernabe B, Weiss MS, Yarrington ME, Barbolina MV, 
# Pannier AK, Jeruss JS, Broadbelt LJ, Shea LD.
# Biotechnol Bioeng. 2011 Feb;108(2):395-403.

##################################################################################################################

# Authors: Beatriz Penalver Bernabe/Dennis Bluver

# Date of Creation: 05/20/2012

# Last modification: 10/21/2013

##################################################################################################################

# NOTE
	
# You will need to download the following packages:
#	-gplots (go to Packages/Install Package(s))
#	-limma (from Bioconductor,www.bioconductor.org)
#	-betr (from Bioconductor,www.bioconductor.org)
#	-sciplot(go to Packages/Install Package(s))

#################################################################################################################

# CLEAR ALL FILES AND FOLDERS AND SET GRAPHICS DEVICE

##Do not adjust

rm(list=ls(all=TRUE))

#################################################################################################################
###					INPUTS												#####
#################################################################################################################

##WORKING DIRECTORY
#################################################################################################################

#Point this to the directory containing all of the data. The "~" represents
#your home folder on a Mac.

	directory = "~/coderuns"

#Name of the folder where you want to save your files

	folder.output = "Stiffness_FINAL"

# Number of factors
 
	nFactor=1

# Factor names. 
# 	For example, if your factors are cell type and drug type, you can name then in the
#	following fashion namesFactors<-c("Cell", "Drug"). Please do not use "treatment"

	namesFactors<-c("Stiffness")

# Number of levels. 
#	For instance, you have 2 cell types and 4 drugs it would be levFactor<-c(2,4)

	levFactor<-c(3)
	
# Colours in which the levels should appear in graphs. If you want to leave by default, just indicate NULL.
# Note it will be in alphabetical order

	valueColors=c("red","green","blue")

# Number of transcription factors 
#	Including "Positive Control", "Normalization Control" (or TA) and "No DNA" 
# 	as separate transcription factors

	nTF = 58

# Normalization vector (TA for instance)

	tfControl<-"TA"

# No DNA or no infected cells

	tfBlank<-"blank"

# Timepoints that you imagine the arrays (e.g. 0,1,2,3)

	timeVector <- c(3,6,9,12,27)

# Maximum number of repeats, number of wells

	nrepeat= 12

# Usual number of repeats-the number of repeats that you have been using for the majority
# of your TFs
 
	nrepeatave=3

# Number of plates
#	Are your TFs split into more than 1 plate that were set up at the same time? If you have more TFs
#	that you could place into one plate, please indicate how many you did use.

	nplate=1

# Number of experiments
#	How many times have you repeated the experiment? If you have different TFs per plate (TFdis="DIF"), how many
#	different days have you plated your arrays?

	nexp=23

# TF distrubution
#	Does all the plates have the same TFs? If yes, please enter "SAME". If you have different 
#     transcription factors in each plate or combination of plates, please enter "DIF"

	TFdis="DIF"

# TF to remove
#	Is there any TF that you would like to remove? Imagine you would like to remove STAT1, then you would indicate
#	TF2remove<-"STAT1". If 2 or more, it would be tf2remove<-c("STAT1", "CRE"). If you do not wish to remove any
# 	indicate TF2remove<-NULL

	TF2remove<-NULL
	
# Specific TF experiment to remove
#	Is there any TF that is not behaving adequately in an specific experiment you would like
# to remove? Imagine you would like to remove STAT1 from exp 3 and CRE from 4 and 
# CRE from 7.
# Here it is how you would indicate. You should create a text file and write two colums. First column
# indicates the TF and the second the experiment
#	STAT1	3
#	STAT1	5
#	CRE		4
#	YY1		1	
# If you do not wish to remove any please indicate TFExp2remove<-NULL

	TFExp2remove<-TRUE
	
# Minimum Experiments a TF should be present to be further consider

	minexp<-3

# Experiments to remove
#	If after see the results from dataAllFluc/dataAllGluc, you decide that some other experiment should be removed
#	please indicate it here (indicate them in the order of reading, i.e, exp2remove<-c(1,5))

	exp2remove<-c(9,10)

# Block per experiment
# When performing the statistical analysis, you have the alternative to block per experiment. This is highly 
# recommended when you do not have a true time 0 or when your experiments do not repeat very well-your correlation
# in the pairs plot is very low (less than 0.7). 

	blockexp<-FALSE

##FILES NAMES
#################################################################################################################

##INSTRUCTIONS FOR SETTING FILE NAMES:

#	If files are directly inside of the working directory (defined above), then just include the
#	names of the files. If the files are within another folder within the working directory,
#	include the full path of the file relative to the working directory.

#	File names for Fluc and Gluc must be set in the following format
#	Experiment1_Time1_plate1, Experiment1_Time1_plate2,...,Experiment1_Time2_plate1,Experiment1_Time2_plate2,...,
#	Experiment1_Time3_plate1,Experiment1_Time3_plate2,...,Experiment2_Time1_plate1,Experiment2_Time1_plate2,...,
#	Experiment2_Time2_plate1, Experiment2_Time2_plate2,...

# Filenames for Gluc (If there are no Gluc files, just type NA)

	filesGluc <-NA

# Filenames for Fluc

	filesFluc<- c("8_31_11/3h_083111.txt","8_31_11/6h_083111.txt","8_31_11/9h_083111.txt","8_31_11/12h_083111.txt","8_31_11/27h_083111.txt",
				"9_6_11/3h_090611.txt","9_6_11/6h_090611.txt","9_6_11/9h_090611.txt","9_6_11/12h_090611.txt","9_6_11/27h_090611.txt",
				"9_28_11/3h_092811.txt","9_28_11/6h_092811.txt","9_28_11/9h_092811.txt","9_28_11/12h_092811.txt","9_28_11/27h_092811.txt",
			     "10_11_11/3h_101111.txt","10_11_11/6h_101111.txt","10_11_11/9h_101111.txt","10_11_11/12h_101111.txt","10_11_11/27h_101111.txt",
			     "10_19_11/3h_101911.txt","10_19_11/6h_101911.txt","10_19_11/9h_101911.txt","10_19_11/12h_101911.txt","10_19_11/27h_101911.txt",
				"10_26_11/3h_102611.txt","10_26_11/6h_102611.txt","10_26_11/9h_102611.txt","10_26_11/12h_102611.txt","10_26_11/27h_102611.txt",
				"11_1_11/3h_110111.txt","11_1_11/6h_110111.txt","11_1_11/9h_110111.txt","11_1_11/12h_110111.txt","11_1_11/27h_110111.txt",
			     "11_9_11/3h_110911.txt","11_9_11/6h_110911.txt","11_9_11/9h_110911.txt","11_9_11/12h_110911.txt","11_9_11/27h_110911.txt",
				"11_21_11/3h_112111.txt","11_21_11/6h_112111.txt","11_21_11/9h_112111.txt","11_21_11/12h_112111.txt","11_21_11/27h_112111.txt",
				"12_1_11/3h_120111.txt","12_1_11/6h_120111.txt","12_1_11/9h_120111.txt","12_1_11/12h_120111.txt","12_1_11/27h_120111.txt",				
				"12_8_11/3h_120811.txt","12_8_11/6h_120811.txt","12_8_11/9h_120811.txt","12_8_11/12h_120811.txt","12_8_11/27h_120811.txt",
				"12_20_11/3h_122011.txt","12_20_11/6h_122011.txt","12_20_11/9h_122011.txt","12_20_11/12h_122011.txt","12_20_11/27h_122011.txt",
				"1_27_12/3h_012712.txt","1_27_12/6h_012712.txt","1_27_12/9h_012712.txt","1_27_12/12h_012712.txt","1_27_12/27h_012712.txt",
				"2_10_12/3h_021012.txt","2_10_12/6h_021012.txt","2_10_12/9h_021012.txt","2_10_12/12h_021012.txt","2_10_12/27h_021012.txt",
				"2_17_12/3h_021712.txt","2_17_12/6h_021712.txt","2_17_12/9h_021712.txt","2_17_12/12h_021712.txt","2_17_12/27h_021712.txt",
				"9_4_12/3h_090412.txt","9_4_12/6h_090412.txt","9_4_12/9h_090412.txt","9_4_12/12h_090412.txt","9_4_12/27h_090412.txt",
				"9_24_12/3h_092412.txt","9_24_12/6h_092412.txt","9_24_12/9h_092412.txt","9_24_12/12h_092412.txt","9_24_12/27h_092412.txt",
				"11_12_12/3h_111212.txt","11_12_12/6h_111212.txt","11_12_12/9h_111212.txt","11_12_12/12h_111212.txt","11_12_12/27h_111212.txt",
				"11_16_12/3h_111612.txt","11_16_12/6h_111612.txt","11_16_12/9h_111612.txt","11_16_12/12h_111612.txt","11_16_12/27h_111612.txt",
				"1_21_13/3h_012113.txt","1_21_13/6h_012113.txt","1_21_13/9h_012113.txt","1_21_13/12h_012113.txt","1_21_13/27h_012113.txt",
				"4_29_13/3h_042913.txt","4_29_13/6h_042913.txt","4_29_13/9h_042913.txt","4_29_13/12h_042913.txt","4_29_13/27h_042913.txt",
				"5_21_13/3h_052013_A.txt","5_21_13/6h_052013_A.txt","5_21_13/9h_052013_A.txt","5_21_13/12h_052013_A.txt","5_21_13/27h_052013_A.txt",
				"07_26_13/072613_3h.txt","07_26_13/072613_6h.txt","07_26_13/072613_9h.txt","07_26_13/072613_12h.txt","07_26_13/072613_27h.txt")

# Condition Maps (include on per plate you have analyzed)

	condMap <-c("8_31_11/condmap_083111.txt","9_6_11/condmap_090611.txt","9_28_11/condmap_092811.txt","10_11_11/condmap_101111.txt",
			"10_19_11/condmap_101911.txt","10_26_11/condmap_102611.txt","11_1_11/condmap_110111.txt","11_9_11/condmap_110911.txt",
			"11_21_11/mech_condmap_112111.txt","12_1_11/mech_condmap_120111.txt","12_8_11/mech_condmap_120811.txt",
			"12_20_11/mech_condmap_122011.txt","1_27_12/condmap_012712.txt","2_10_12/mech_condmap_021012.txt",
			"2_17_12/mech_condmap_021712.txt","9_4_12/condmap2_090412.txt","9_24_12/condmap_092412.txt","11_12_12/conditionmap_111212.txt",
			"11_16_12/conditionmap_111612.txt","1_21_13/condmap_012113.txt","4_29_13/condmap_042913.txt","5_21_13/condmapA_052013.txt",
			"07_26_13/condmap_072613.txt")

# Condition Lists (include on per plate you have analyzed)
#
	condList <- c("8_31_11/condlist_083111.txt","9_6_11/condlist_090611.txt","9_28_11/condlist_092811.txt","10_11_11/condlist_101111.txt",
			"10_19_11/condlist_101911.txt","10_26_11/condlist_102611.txt","11_1_11/condlist_110111.txt","11_9_11/condlist_110911.txt",
			"11_21_11/mech_condlist_112111.txt","12_1_11/mech_condlist_120111.txt","12_8_11/mech_condlist_120811.txt",
			"12_20_11/mech_condlist_122011.txt","1_27_12/MECHcondlist_012712.txt","2_10_12/mech_condlist_021012.txt",
			"2_17_12/mech_condlist_021712.txt","9_4_12/mech_condlist_090412.txt","9_24_12/MECH_condlist_092412.txt",
			"11_12_12/conditionlistMECH_111212.txt","11_16_12/conditionlistMECH_111612.txt","1_21_13/conditionlistMECH_012113.txt",
			"4_29_13/mech_condlist_042913.txt","5_21_13/MECHcondlistA_052013.txt","07_26_13/condlist_MECH_072613.txt")

##NORMALIZATION AND STATISTICAL PARAMETERS
#################################################################################################################

#Removal of data points below the background
#	If you want to eliminate all the datapoints below the background, please indicate bg=1
	
	bg=1

#Confidence interval for background data

     confBg=0.995

#Do you have a positive control? If yes, indicate the name. Otherwise indicate tf.pc<-NA
#Positive control name

	tfPc<-NA

#Confidence interval for positive control above the TA

     confPc=0.995

# Number of times about the background for positive control

	aboveBcPc<-3

#Removal of outliers
#	If you want to eliminate all the outlier, please indicate outlier.removal=1

	outlierRemoval=1

#Confidence interval for outlier removal (we recommend a high number so that you are not removing 
# data that are not outlier due to the skewness of the data or the limit amount of repeats available

     confOutlier=0.9995

# Number of experimental repeats for all the transcription factors. If you are uncertain, write NULL 

	expTF<-NULL

# Scale free the data.If yes, indicate TRUE, otherwise indicate FALSE

	scaleFree=TRUE
	
# Scale the mean (New after version 5.0). This option will just substract the average mean of the entire
# experiment, so that all the experiments for a given TFs would have the same absolute mean. This is
# mostly helpful for NLCA.

	scaleMean=FALSE

# Normalize by the initial timepoint-in case that all the first time points are supposed to be identical
# (i.e, if all the cells were measured at time 0 before giving the treatment)If yes, indicate TRUE, otherwise
# indicate FALSE

	t0Norm=FALSE

# If you have gluc data available, would you like to use the initial Gluc value for each treatment to normalize 
# the Fluc data? It can help as an indication of the total number of cells in each well. If yes, indicate TRUE, otherwise
# indicate FALSE

	tGluc=FALSE
	
# Would you like to remove some of the experimental variability by multiplying e
# each experimental curve by the median of the initial value?

	tWithin=FALSE

# Do you want to shift the curves to the initial timepoint? It would assume that all the experimental repeats
# should have the same value for a given treatment-remove experimental variability. If you do not want to use 
# this option please indicate tin=0. If you would like to remove experimental variability, please indicate the 
# following:
#	if all the initial timepoints should be the same-for instance, you measure at time 0 before adding the
#		different treatments, please indicate tin=1. 
#	if all the initial timepoints are not the same for the different treatments-for instance, your initial
#		timepoints start after 1 day of treatment or you are doing experiments in 2D and 3D, please indicate 
#		tin=2
#	if you have more than 1 factor, please indicate in this format, tin<-c(1,2). In this case, you would like
#		to shift both factors, where the first factor all the treatment should be the same value for the initial
#		timepoint but not for the second.
	
	tIn<-0
	
# Untreated condition or reference you would like to use for comparison. If none specifically, indicate 
# nameUntreated<-"none". If you have several factors, for instance, cells type MCF7 and drug NT, please separate
# them by , and quoted, i.e, nameUntreated<-c("MCF7","NT")

	nameUntreated<-c("none")

# Custom contrasts- it is not currently implemented
 
	customcontrast<-NULL

# Number of bootstrapping iterations if required

	nboot=200000

# Method for fdr correction. You can select "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
# "fdr" or "none"

  	methodlimma<-"fdr"
 
# p-value cut off

	pvaluelimma<-0.001
	
# Use meta-analysis

	metanalysis<-TRUE

# p-value cut-off metanalysis

	pvaluemetanalysis<-0.001

##OTHER PARAMETERS
#################################################################################################################

# Indicate your Operating Systems
# 	"PC" or "Mac" 
	
	OS<-"Mac"

# Location of the source

	sourcelocation<-"~/coderuns/functionsALCAv5.1.R"
																		  
	
#################################################################################################################
####																		 ####
	####																######
		####														####
	####																####
###																		####
#################################################################################################################


#################################################################################################################
###					READING AND ORGANIZING INPUT DATA								#####
#################################################################################################################
#################################################################################################################

#	Each file is read, the outliers are removed and it is normalized by TA


# INITIALIZE VARIABLES
##############################################################################################################

inputParameters<-list(
    directory=directory,
 	nFactor=nFactor, 
	namesFactors=namesFactors,
	levFactor=levFactor,
	valueColors=valueColors,
	nTF = nTF,
	tfControl=tfControl,
	tfBlank=tfBlank,
	timeVector=timeVector,
	nrepeat= nrepeat,
	nrepeatave=nrepeatave,
	nplate=nplate,
	TFdis=TFdis,
    tfloc="FIX",
	nexp=nexp,
	TF2remove=TF2remove,
	TFExp2remove=TFExp2remove,
	filesGluc=filesGluc,
	filesFluc=filesFluc,
	condMap=condMap,
	condList=condList,
	bg=bg,
	confBg=confBg,
    tfPc=tfPc,
	confPc=confPc,
	aboveBcPc=aboveBcPc,
	outlierRemoval=outlierRemoval,
	confOutlier=confOutlier,
	expTF=expTF,
	scaleFree=scaleFree,
	scaleMean=scaleMean,
	t0Norm=t0Norm,
	tGluc=tGluc,
	tWithin=tWithin,
	tIn=tIn,
	nameUntreated=nameUntreated,
	OS=OS,
	sourcelocation=sourcelocation,
	exp2remove=exp2remove,
	minexp=minexp,
	blockexp=blockexp,
	nboot=nboot,
	customcontrast=customcontrast,
	methodlimma=methodlimma,
	pvaluelimma=pvaluelimma,
	folder.output=folder.output,
	metanalysis=metanalysis,
	pvaluemetanalysis=pvaluemetanalysis)

# Open source

	source(inputParameters$sourcelocation)	
	
# Set working directory

	setwd(inputParameters$directory)

# Create directory

     dir.create(folder.output)

# INITIALIZE TABLES AND VECTORS
##############################################################################################################
  
# Initialize time, experiments and plate counter

  time=0
  exp=0
  plate=0

# Determine factor, levels and treatments

  UpdatedInput<-factorlist(inputparameters=inputParameters)
  inputParameters<-UpdatedInput$inputparameters
  expCond<-UpdatedInput$ExpCond

# READ FILES
##############################################################################################################
  if (!is.na(filesGluc[1])){
	dataAllGluc<-tableGen(inputparameters=inputParameters,expCond=expCond,comp="Gluc")
  } else {
	dataAllGluc<-NA
  }
  #Generation of the tables, background and outlier correction
  dataAllFluc<-tableGen(inputparameter=inputParameters,expCond=expCond,comp="Fluc")
  
  # Plot raw data in linear and log scale
  plotRawDataPDF(dataAll=dataAllFluc,inputparameters=inputParameters,
	expCond=expCond,comp="Fluc")
  
  # Plot raw data after background and outlier removal in linear and log scale
	plotpreNormData(dataAll=dataAllFluc,inputparameters=inputParameters,
	expCond=expCond,comp="Fluc")

# NORMALIZATION
###############################################################################################################################################################

   #Normalize the data by the initial average mean per repeat (optional) and
   #normalize by TA
   dataAllFluc<-withInAndTA(dataGluc=dataAllGluc,dataFluc=dataAllFluc,
	inputparameters=inputParameters,expCond=expCond)
	
	# Normalize the data by all the selected options, such as Gluc, t0,
	# shift the data,...
    dataAllFluc<-normBetween(dataFluc=dataAllFluc,
	inputparameters=inputParameters,expCond=expCond)
	
	# Plot all the normalization figures
    plotNormalization(dataAll=dataAllFluc,inputparameters=inputParameters,
		expCond=expCond)

   PlotBoxPlot(dataAll=dataAllFluc,inputparameters=inputParameters,
	expCond=expCond,cex.main=0.8,cex.lab=0.8,cex.axis=0.8)

  if (nexp>1){
  	PlotPairs(dataAll=dataAllFluc,
		inputparameters=inputParameters,expCond=expCond,
		cex.main=0.8,cex.lab=0.8,cex.axis=0.8)
  }

# STATISTICAL ANALYSIS
###############################################################################################################################################################

  dataAllFluc<-tableOrganization(dataFluc=dataAllFluc,
	inputparameters=inputParameters,expCond=expCond)

  dataAllFluc<-statAnalysis(dataFluc=dataAllFluc,
	inputparameters=inputParameters,expCond=expCond)

