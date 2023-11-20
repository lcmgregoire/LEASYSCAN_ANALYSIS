##################
# EZTr_main #

# This framework is developed as part of the research article:
# Automated Discretization of ‘Transpiration Restriction to Increasing VPD’ Features from Outdoors High-Throughput Phenotyping Data
#
# Soumyashree Kar, Ryokei Tanaka, Lijalem Balcha Korbu, Jana Kholová, Hiroyoshi Iwata, Surya S. Durbha, J. Adinarayana, Vincent Vadez
##################

#' HTP data processing
#'
#' Convert  gravimetric sensors or load cells (LC) time series into
#' Evapotranspiration (ETr) and Transpiration (Tr) time series.
#'
#' The function iterates over the different time points (at 15 min interval) of
#' each sector or pot to generate a matrix of Tr values. Fifteen biologically
#' relevant features are extracted from the Tr time series for each day,
#' followed by computation of feature heritability.
#'
#' The entire process comprises four steps: 1st - Conversion of LC data to ETr;
#' 2nd - Calculation of reference i.e. Penman Monteith ET for the given weather
#' condition, and filtering ETr based on reference ET; 3rd - Generating smooth
#' ETr time series; 4th - Extraction of Tr from smooth ETr, Tr features and
#' each feature's broad-sense heritability estimate.
#' 
#' Input Files required:
#' 1.Loadcells
#' 2.Loadcells metadata
#' 3.Sensor climate data
#' 4.Sensor unit map
#' 5.Plant eye
#'
#' #' Functions used for each step:
#' Step1 - extractRawLCmatrix(), curateRawLCgetET(), filterETrExtremeCols();
#' Step2 - extractWthrVar(), prepcsWthr(), calculateETref(),generateETrRatio(),
#'         genThreshVal(), dataPART(), threshETr(), ordFiltETr();
#' Step3 - smoothETr();
#' Step4 - calculateTr(), getFeatures(), getFeatureHe().
#'
#'
#' The function includes seven inputs, 'data' includes LC data,
#' experimental design and genotype-replicate information in metadata, weather
#' data and genotype-replicate specific leaf area data; first (Date1) and
#' last (Date2) dates of the experiment; 'irrg.dts' is a vector of dates when
#' plants were irrigated; 'opPATH' and opPATH.smth to store intermediate
#' feature-specific files.
#'
#' @param HTP_data \code{list} of 5 dataframes composed of input files
#' a) loadcell data with unit, genotype, g_alias, treatment, timestamp
#' and Mass(g) columns
#' b) metadata with unit, old_unit, Experiment, Treatment, Species
#' Genotype, G. Alias and Replicates columns
#' c) weather data with sensor, variable, timestamp, value columns
#' d) solar radiation data is stored in separate dataframe in format c)
#' e) leaf area data with Sector, Experiment, Treatment, Species,
#' Genotype, G..Alias, Replicates, timestamp and leaf area columns
#'
#' @param lastDate \code{string} containing last date of experiment
#' in YYYY-MM-DD format.
#'
#' @param irrg.dts \code{vector} of irrigation dates, each date mentioned as
#' string in YYYY-MM-DD format. If no such date is needed, same as the last date.
#'
#' @param Date1 \code{string} value specifying first date of time series in
#' YYYY-MM-DD hh:mm:ss format.
#'
#' @param Date2 \code{string} value specifying last date of time series in
#' YYYY-MM-DD hh:mm:ss format.
#'
#' @param opPATH \code{string} value specifying the directory path for
#' intermediate results to be stored.
#'
#' @param opPATH.smth \code{string} value specifying the directory path for
#' feature-specific results to be stored.
#'
#'
#' @return Return:
#'
#' The function returns a list of fifteen outcomes:
#' LCraw_tsmeta = Matrix of time stamp values present in raw load cell data,
#' LCraw = Matrix of raw load cell data,
#' ETrmm_Obs = Matrix of raw or observed ETr in mm
#' ETr_Obs_ERR.SEC.nms = Matrix of sector names with very high proportion of extreme values,
#' ETmm_Obs_FINAL = Matrix of ETr in mm after removing the erroneous sectors,
#' Wthr_agg15min = Matrix of all weather variables aggregated for 15 minutes interval,
#' Wthr.ETref.ETobs = Matrix of combined time series of Reference ET, ETr, weather values,
#' ETrRatio_TW1_ThreshVALS = Matrix of filtered values of Day-time ref ET/ETr ratios,
#' ETrRatio_TW2_ThreshVALS = Matrix of filtered values of Night-time ref ET/ETr ratios,
#' IrrgFilt_ETr = Matrix of Irrigation filtered ETr matrix,
#' IrrgFilt_ETr_Imptd = Matrix of irrigation filtered ETr matrix imputed,
#' ETr_smth = Matrix of smooth ETr time series,
#' Tr = Matrix of Tr time series,
#' featureH2 = Matrix of each feature heritability estimate on each day,
#' eachFeature_TS = List of each feature's time series for all genotypes.
#'
#' @author Soumyashree Kar email<ksoumya2301@gmail.com>
#'
#' @references
#'
#' Vadez, V., Kholová, J., Hummel, G., Zhokhavets, U., Gupta, S. K., &
#' Hash, C. T. (2015). LeasyScan: a novel concept combining 3D imaging and
#' lysimetry for high-throughput phenotyping of traits controlling plant water
#' budget. Journal of Experimental Botany, 66(18), 5581-5593.
#' 
#' 
#install.packages("easypackages")

setwd(dir = "C:/Users/2021lg003/Documents/Sorghum/Experiments/EXP_INDE_2023/ANALYSIS_ALL_DATASET/DATA_SMOOTHING")

#install.packages("readxl", "hms", "xts", "dplyr","mgcv","PerformanceAnalytics", "wavelets",  
               #  "signal", "tidyverse", "zoo", "h2o", "sqldf", "ggplot2", "plyr",
                # "lubridate", "BioFTF", "plantecophys", "highfrequency", "stringr",
                # "chron", "nonlinearTseries", "tsfeatures", "splitstackshape", "psych")
require(devtools)
install_version("lubridate", version = "1.7.1", repos = "http://cran.us.r-project.org")
library(easypackages)
#update.packages("easypackages")

library("readxl")
library("hms")
library("xts")
library("dplyr")
library("mgcv")
library("PerformanceAnalytics")
library( "wavelets")
library("signal")
library( "tidyverse")
library( "zoo")
library( "h2o")
library( "sqldf")
library("ggplot2")
library("plyr")
library("BioFTF")
library("plantecophys")
library( "highfrequency")
library( "stringr")
library("chron")
library( "nonlinearTseries")
library("tsfeatures")
library("splitstackshape")
library("psych")

packageVersion("lubridate") 

library(lubridate)
# Load all the functions #
# Functions needed for processing Stage-I: LC to ETr extraction
source('./functions/extractRawLCmatrix.R') # Included 'Treatment col in meta.d'
source('./functions/curateRawLC.R')
source('./functions/filterLCExtremeCols.R')
source('./functions/getETr.R')
source('./functions/convETr.R')
source('./functions/filterETrExtremeCols.R')

# Functions needed for processing Stage-II: ETref extraction and ETr thresholding
source('./functions/extractWthrVar.R')
source('./functions/prepcsWthr.R')
source('./functions/calculateETref.R')
source('./functions/generateETrRatio.R')
source('./functions/genThreshVal.R')
source('./functions/dataPART.R')
source('./functions/threshETr.R')
source('./functions/ordFiltETr.R')

# Functions needed for processing Stage-III: raw and smooth Tr, Tr-features and feature-H2 extraction
source('./functions/calculateTr.R')
source('./functions/getFeatures.R')
source('./functions/getFeatureHe.R')
source('./functions/getFeatures(1).R')
source('./functions/smoothETr.R')
source('./functions/getFeatures.R')

# Load data
load("data/PERIODE2-01-10-10-10/IRD_data_Exp_NEW.RData")
allData <- IRD_data_Exp_NEW


ldt <- readline(prompt = "Enter LAST DATE of experiment (YYYY-MM-DD): ")
lastDate=as.Date(as.character(ldt)) #2023-10-10

fdt <- readline(prompt = "Enter FIRST DATE of experiment (YYYY-MM-DD): ")
firstDate=as.Date(as.character(fdt)) #2023-10-01

seq_by <- readline(prompt = "Enter desired time (in min) interval, e.g. 15/30/45/60: ")

# NOTE: Enter date(s) of irrigation or when data is extremely noisy.
# NOTE: If no such date, enter LAST DATE of experiment
irrg.dts <- c("2023-10-10")

Date1="2023-09-30 23:45:00" # Day before 'firstDate'(For 15min, "2021-04-18 23:46:00")
Date2="2023-10-10 23:45:00" #'lastDate'(For 15min, "2022-04-29 23:45:00")

opPATH <- "./results/"    #####opPATH <-"C:/Users/2021lg003/Documents/Manip_loadcells/Lissage_Manip_BenoitClerget_092022/results/"
opPATH.smth="./results/smthFeaturesTimeSeries/"   ####opPATH.smth="C:/Users/2021lg003/Documents/Manip_loadcells/Lissage_Manip_BenoitClerget_092022/results/smthFeaturesTimeSeries"
opPATH.raw="./results/rawFeaturesTimeSeries/" ###opPATH.raw="C:/Users/2021lg003/Documents/Manip_loadcells/Lissage_Manip_BenoitClerget_092022/results/rawFeaturesTimeSeries/"
opPATH.smthTR="./results/smthTR_FeaturesTimeSeries/"  ###opPATH.smthTR="C:/Users/2021lg003/Documents/Manip_loadcells/Lissage_Manip_BenoitClerget_092022/results/smthTR_FeaturesTimeSeries/" 


# Get load cells data i.e. weights of sector
m.lc <- allData$m.lc

# Get Genotype and Exp design metadata
meta.d <- allData$meta.d
meta.d <- distinct(meta.d) #remove duplicate entries
meta.d <- meta.d[,1:8]

# Get the list of species from the metadata
species.nm <- unique(meta.d$Species)


### IMP- For mulch dataset repeat below for each species ###

# Include the species ID e.g. 1 for 'Oryza' as per the data
meta.d.sp <- meta.d[meta.d$Species==species.nm[1], ]

# Find sectors with missing metadata
noEntrySecs <- which(!unique(m.lc$unit) %in% unique(meta.d.sp$unit))

noEntrySecNms <- unique(m.lc$unit)[noEntrySecs]

# Remove sectors-without-metadata from original loadcells file
m.lc <- m.lc[! m.lc$unit %in% noEntrySecNms, ]


####### Stage-I: Process LC data and generate ETr matrix #######
st.time <- Sys.time()

# Extract matrix of loadcell data #

LC.MAT.OP <- extractRawLCmatrix(x = m.lc,   y = meta.d.sp,
                                s=firstDate, z = lastDate, 
                                inter = seq_by)

LC.MAT.f <- LC.MAT.OP$LC.MAT.f

LC.MAT.TSinfo <- LC.MAT.OP$LC_tsmeta


# Add metadata to LC matrix - For Mulch data, add 'Trtmt' #
# select from Metadata: "unit", "old_unit", 'Trtmt', "Genotype, "G..Alias", "Replicates"
meta.d.LCmat <- meta.d.sp[meta.d.sp$unit %in% colnames(LC.MAT.f)[-1], 
                          c(1,2,4,6,7,8)] 

LC.MAT.f.t <- as.data.frame(t(LC.MAT.f))

colnames(LC.MAT.f.t) <- LC.MAT.f$TS

LC.MAT.f.t <- LC.MAT.f.t[-1,]

# reorder rows of 'meta.d.sp' according to rownames/unit of LC.MAT.f.t
meta.LCDF <- meta.d.LCmat[order(match(meta.d.LCmat$unit, rownames(LC.MAT.f.t))), ]

LC.MAT.raw <- as.data.frame(cbind(meta.LCDF, LC.MAT.f.t))

write.csv(LC.MAT.raw, paste0(opPATH, "OP-1","_LCraw_wNA.csv"))


# Start outlier detection, removal and imputation of LC Matric to generate ETr profiles #
# Pre-process raw LC data: outliers removal and imputation #
imputed.DF.final <- curateRawLC(x = LC.MAT.f, y = meta.LCDF)

write.csv(imputed.DF.final, paste0(opPATH,"OP-2","_LC_olrm_imputed.csv"))


# Identify the highly extreme valued sectors #
err.sec.info <- filterLCExtremeCols(x = imputed.DF.final, y = meta.LCDF)

err.sec.nm <- err.sec.info$err.sec.NM

err.sec.meta <- err.sec.info$err.sec.META

write.csv(err.sec.meta, paste0(opPATH, "OP-3","_LCimp_errorUnits.csv"))


# Remove the err.cols i.e. sectors with extreme values
impData.errSEC.rmvd <- imputed.DF.final[!imputed.DF.final$unit %in% err.sec.nm, ]

write.csv(impData.errSEC.rmvd, paste0(opPATH,"OP-4","_LCimp_ETr_IP.csv"))


# Generate ETr profiles from "impData.errSEC.rmvd" dataframe #
et.vals <- getETr(x = impData.errSEC.rmvd)

et.obs <- et.vals$obsETr_core

ETr_Meta <- et.vals$obsETr_meta


# Convert ETr in grams to mm (Y/N) #

ETr_F <- convETr(x = ETr_Meta, y = et.obs)
write.csv(ETr_F, paste0(opPATH, "OP-5","_ETr_Obs.csv"))


# Identify error plots from ETr values using the similar method as above #
ETr_err.sec.info <- filterETrExtremeCols(x = ETr_Meta, y = meta.LCDF)

err.sec.nm <- ETr_err.sec.info$ETr_err.sec.NM

err.sec.meta <- ETr_err.sec.info$ETr_err.sec.META

if(length(err.sec.nm)>0){write.csv(err.sec.meta,  
                                   paste0(opPATH, "OP-6","_ET_Obs_ERR.SEC.nms.csv"))}


# Remove the err.cols i.e. sectors with extreme values #
ETr_Meta_ERRsec.rmvd <- ETr_Meta

ETr_Meta_ERRsec.rmvd <- ETr_Meta_ERRsec.rmvd[!ETr_Meta_ERRsec.rmvd$unit %in% err.sec.nm, ]

write.csv(ETr_Meta_ERRsec.rmvd, paste0(opPATH, 
                                       "OP-7","_ETr_Obs_FINAL.csv"))
##### END ######


####### Stage-II: Process Weather data to obtain ETref and ETr ratio matrix #######

wthr.DFagg15min <- IRD_data_Exp_NEW$climate
wthr.DFagg15min$Date <- lubridate::ymd(wthr.DFagg15min$Date) #Format here

wthr.DFagg15min <- wthr.DFagg15min[wthr.DFagg15min$Date>=firstDate & 
                                     wthr.DFagg15min$Date<=lastDate,]

#Ensure to retain only till Row No.with complete data 
rn <- which(is.na(wthr.DFagg15min$Date)==TRUE)
wthr.DFagg15min <- wthr.DFagg15min[-c(rn[1]:nrow(wthr.DFagg15min)),] 
wthr.DFagg15min$TS <- paste0(as.character(wthr.DFagg15min$Date), " ",
                                    wthr.DFagg15min$TS)
wthr.DFagg15min$TS <- ymd_hms(wthr.DFagg15min$TS) ##pb with lubridate? 
wthr.DFagg15min <- wthr.DFagg15min[,-1]

# Create empty base matrix for weather variables
end.seq_hhmm  <- '23:45'

if (seq_by == '15'){
  end.seq_hhmm  <- '23:45'
} else if(seq_by == '30'){
  end.seq_hhmm  <- '23:30'
} else if(seq_by == '45'){
  end.seq_hhmm  <- '23:15'
} else if(seq_by == '60'){
  end.seq_hhmm  <- '23:00'}

TS_base<-as.data.frame(as.character(seq(ymd_hm(paste0(firstDate," ",'00:00')),
                                        ymd_hm(paste0(lastDate," ",end.seq_hhmm)), 
                                        by = paste0(seq_by, " ", "mins"))))
#  
# TS_base<-as.data.frame(as.character(seq(ymd_hm(paste0(firstDate," ",'00:00')),
#                                         ymd_hm(paste0(lastDate," ",'23:45')), by = '15 mins')))

names(TS_base)[1]<-c("int.val")

TS_base$time <- strftime(TS_base$int.val, format="%H:%M:%S", tz="UTC")

# Since, the hms values slightly differ in the orginal dataset than the ideal 15min interval values,
# replace them with the TS_base strftime (hms) values

hms.ts.base<-unique(TS_base$time)
names(hms.ts.base)<-c("time")

print("Climate matrix timestamp mapping status")

i<-nrow(wthr.DFagg15min)
pbar <- create_progress_bar('text')
pbar$init(i)

for(i in 1:nrow(wthr.DFagg15min))
{
  if(! wthr.DFagg15min$TS[i] %in% hms.ts.base)
  {
    j <- which.min(abs(chron(times=hms.ts.base) - 
                         chron::chron(times=format(wthr.DFagg15min$TS[i],format = "%H:%M:%S")))) # find which value in ts-base-vector is nearest to each-DP
    
    # assign the nearest ts-base-vector value to that DP
    dd <- date(wthr.DFagg15min$TS[i])
    wthr.DFagg15min$TS[i] <- ymd_hms(paste0(dd, TS_base$time[j]))} 
  pbar$step()
}

# TS_ALL<-as.data.frame(as.character(seq(ymd_hm(paste0(firstDate," ",'00:00')),
#                                        ymd_hm(paste0(lastDate," ",'23:45')), by = '15 mins')))
TS_ALL<-as.data.frame(as.character(seq(ymd_hm(paste0(firstDate," ",'00:00')),
                                       ymd_hm(paste0(lastDate," ",end.seq_hhmm)), 
                                       by = paste0(seq_by, " mins"))))
names(TS_ALL)[1]<-c("TS.n") # MUST be the same as in original data set m.lc.df

TS_ALL$TS.n <- ymd_hms(TS_ALL$TS.n)

# create empty dataframe to store all values

Wthr.MAT <- as.data.frame(matrix(nrow = length(TS_ALL$TS.n), 
                               ncol = ncol(wthr.DFagg15min)))

Wthr.MAT[ ,1] <- TS_ALL$TS.n; names(Wthr.MAT)[1] <- "TS"
colnames(Wthr.MAT)[2:ncol(Wthr.MAT)] <- names(wthr.DFagg15min)[2:ncol(wthr.DFagg15min)]

df <- merge(x = wthr.DFagg15min, y = Wthr.MAT, by = "TS", all = TRUE) # perform outer join to merge by id=TS.n
df.new <- df[,1:6]; names(df.new)[2:6]<-names(wthr.DFagg15min)[2:6]
df.new <- df.new[!duplicated(df.new[c("TS")]),]

etrDTS <- as.data.frame(colnames(ETr_Meta_ERRsec.rmvd)[7:ncol(ETr_Meta_ERRsec.rmvd)])
names(etrDTS) <- "TS"
etrDTS$TS <- ymd_hms(etrDTS$TS)

df.new <- df.new[df.new$TS %in% etrDTS$TS, ] # subset after outer join to ensure that NAs don't add extra rows

# Pre-process weather
df.new$Temp <- if(sum(is.na(df.new$Temp))>0){df.new$Temp<-na.aggregate.default(df.new$Temp)}
df.new$RH <- if(sum(is.na(df.new$RH))>0){df.new$RH<-na.aggregate.default(df.new$RH)}
df.new$SR <- if(sum(is.na(df.new$SR))>0){df.new$SR<-na.aggregate.default(df.new$SR)}
df.new$WS <- if(sum(is.na(df.new$WS))>0){df.new$WS<-na.aggregate.default(df.new$WS)}
wthr.DFagg15min.filt <- df.new ##pb juste au dessus -- enleve les colonnes alors qu'il n'y a pas de NA..? 

# Compute VPD and insert into the weather DF #
SVP <- 610.7*(10^(7.5*wthr.DFagg15min.filt[ ,2]/(237.3+wthr.DFagg15min.filt[ ,2])))
VPD <- ((1 - (wthr.DFagg15min.filt[ ,3]/100))*SVP)/1000
wthr.DFagg15min.filt[ ,4] <- VPD

et.obs <- ETr_Meta_ERRsec.rmvd


###################################################

# Calculate Penman Monteith ET #
wthr.df1 <- calculateETref(x=wthr.DFagg15min.filt)
wthr.ETref.df <- as.data.frame(wthr.df1)
empty.MAT <- matrix(nrow = 8, 
                    ncol = (ncol(et.obs)-nrow(wthr.df1)))

# select columns "Temp"  "RH"    "VPD"   "SR"    "WS"    "Tmax"  "Tmin"  "ETref"
empty.MAT.wthr.ETref <- as.data.frame(cbind(empty.MAT, t(wthr.ETref.df[,c(2:6, 9:11)])))
colnames(empty.MAT.wthr.ETref) <- colnames(et.obs)
wthr.ETref.ETobs <- as.data.frame(rbind(empty.MAT.wthr.ETref, et.obs))

# Remove irrigation dates #
file.colnms <- colnames(wthr.ETref.ETobs)
wthr.ETref.ETobs <- wthr.ETref.ETobs[ ,!substr(file.colnms,1,10) %in% irrg.dts]

write.csv(wthr.ETref.ETobs, paste0(opPATH, "OP-8","_wthr.ETref.ETobs.csv"))


######### Start ET0 Ratio-based Thresholding ########
# ET_ratio_mat <- generateETrRatio(x = wthr.ETref.ETobs)
# 
# write.csv(ET_ratio_mat, paste0(opPATH, "OP-15","_wthr.ETref.ETobs.Ratio.csv"))
# 
# ETr_baseFILE <- wthr.ETref.ETobs
# baseFILE <- ET_ratio_mat
# 
# # Assign "G" to NA Genotype cells
# baseFILE$Genotype[9:nrow(baseFILE)][is.na(baseFILE$Genotype)[9:nrow(baseFILE)]] <- "G" 
# 
# # Make the original date sequence in the LC time series from which irrg. dates will be filtered
# act.dts <- seq(date(wthr.DFagg15min.filt$TS[1]), 
#                date(wthr.DFagg15min.filt$TS[nrow(wthr.DFagg15min.filt)]), 1)
# 
# act.dts.ts <- rep(act.dts, each = 96)
# 
# act.dts.ts.df <- data.frame(org.dts.ts = c(rep(NA, 
#                  (ncol(ETr_baseFILE) - length(act.dts.ts))),
#                                            as.character(act.dts.ts)))
# act.dts.ts.df$org.dts.ts <- as.Date(act.dts.ts.df$org.dts.ts)
# 
# rownames(act.dts.ts.df) <- colnames(ETr_baseFILE)
# 
# # Identify columns from ETr raw data which need to be filtered
# irrg.cols <- which(act.dts.ts.df$org.dts.ts %in% as.Date(irrg.dts))
# 
# # Prepare a copy of raw data as the Input File for the Smoothing process
# ETr_smth_IP <- ETr_baseFILE
# Ratio_smth_IP<- baseFILE
# 
# 
# ETr_smth_IP <- ETr_smth_IP[ ,-irrg.cols]
# Ratio_smth_IP <- Ratio_smth_IP[ ,-irrg.cols]
# 
# 
#         ### Changed group_by "Genotype" to "Treatment" ###
#   ### Check column# accordingly i.e. '3' in part:c(3, 7:ncol(Ratio_smth_IP)) ###
# 
# # Thresholding using 2 time-windows 06:30 to 18:30 and remaining #
# by_Genotype <- Ratio_smth_IP[-c(1:8),c(3, 7:ncol(Ratio_smth_IP))] %>% group_by(Treatment)
# 
# by_Genotype_Mean <- by_Genotype  %>% summarise_all(~mean(.))
# 
# by_Genotype_Mean$Treatment <- factor(by_Genotype_Mean$Treatment)
# 
# geno.ETr <- as.data.frame(t(by_Genotype_Mean))
# 
# geno.ETr <- geno.ETr[-1, ]
# 
# geno.ETr <- as.data.frame(apply(geno.ETr, 2, function(x) {as.numeric(as.character(x))}))
# 
# colnames(geno.ETr) <- paste0("G_", by_Genotype_Mean$Treatment)
# 
# ETref <- as.data.frame(Ratio_smth_IP[8, 7:ncol(Ratio_smth_IP)])
# 
# baseDF <- as.data.frame(cbind(t(ETref), geno.ETr))
# colnames(baseDF)[1] <- c("ETref")
# 
# avg.rawETr <- as.data.frame(cbind(colnames(by_Genotype)[-1],geno.ETr))
# names(avg.rawETr)[1] <- "Timestamp"
# 
# write.csv(avg.rawETr, paste0(opPATH,"OP-16","_rawETr_Layla2a_avg.csv"))
# 
# 
# # Check date format in rownames and run accordingly
# baseDF$TS <- ymd_hms(rownames(baseDF)); 
# baseDF$date <- date(baseDF$TS); 
# baseDF$time <- strftime(baseDF$TS, format="%H:%M:%S", tz="UTC")
# baseDF$solarRAD <- as.numeric(as.character(Ratio_smth_IP[4, 7:ncol(Ratio_smth_IP)])) ## Need to convert row to vector using c()
# 
# baseDF <- baseDF[ ,c((ncol(baseDF)-3), (ncol(baseDF)-2), (ncol(baseDF)-1), ncol(baseDF), 1:(ncol(baseDF)-4))]
# 
# 
# ### Steps to find the Thresholds and filter raw ETr ###
# 
# # 1. First find the set of unique dates.
# # 2. Divide the whole data into 2 time-windows, TW (06:30 - 18:30, Rest).
# # 3. For each TW, get the threshold values.
# # 4. Partition ETr as per TWs and apply thresholds.
# # 5. Filter ETr, impute and merge the TWs.
# 
# 
# # 1. First find the set of unique dates
# unq.dts <- unique(baseDF$date)
# 
# # 2. Divide the whole data i.e. baseDF into 2 time-windows, TW (06:30 - 18:30, Rest)
# base_TW1 <- baseDF[baseDF$time >= "06:30:00" & baseDF$time < "18:30:00", ]
# base_TW2 <- baseDF[!baseDF$time %in% base_TW1$time, ]
# 
# # 3. For each TW, get the threshold values. 
# baseTW1_ThreshVALS <- genThreshVal(x = unq.dts, y = base_TW1)
# baseTW2_ThreshVALS <- genThreshVal(x = unq.dts, y = base_TW2)
# 
# TW1.thresh <- ceiling(median(baseTW1_ThreshVALS$Q_75, na.rm = TRUE))
# TW2.thresh <- ceiling(median(baseTW2_ThreshVALS$Q_75, na.rm = TRUE))
# 
# write.csv(baseTW1_ThreshVALS, paste0(opPATH, "OP-17a","_ETrRatio_TW1_ThreshVALS.csv"))
# write.csv(baseTW2_ThreshVALS, paste0(opPATH, "OP-17b","_ETrRatio_TW2_ThreshVALS.csv"))
# 
# 
# # 4. Make ETr partitions
# ETr_TWs <- dataPART(x = ETr_smth_IP)
# 
# ETr_P1 <- ETr_TWs$P1 
# ETr_P2 <- ETr_TWs$P2
# 
# 
# # Make ETr Ratio partitions
# ETr_Ratio_TWs <- dataPART(x = Ratio_smth_IP)
# 
# ETr_Ratio_P1 <- ETr_Ratio_TWs$P1
# ETr_Ratio_P2 <- ETr_Ratio_TWs$P2
# 
# 
# # Apply TW-specific thresholds
# ETr_filt_TW1<- threshETr(x = TW1.thresh, y = ETr_Ratio_P1, z = ETr_P1)
# ETr_filt_TW2<- threshETr(x = TW2.thresh, y = ETr_Ratio_P2, z = ETr_P2)
# 
# # Save filtered ETr data
# ETr_filt_TW1_ord <- ordFiltETr(x = ETr_filt_TW1, y = ETr_smth_IP)
# ETr_filt_TW2_ord <- ordFiltETr(x = ETr_filt_TW2, y = ETr_smth_IP)
# 
# ETr_filt_ord <- as.data.frame(rbind(ETr_filt_TW1_ord, ETr_filt_TW2_ord))
# 
# ETr_filt_ord$TS <- ymd_hms(rownames(ETr_filt_ord))
# 
# ETr_filt_ord <- ETr_filt_ord[order(ETr_filt_ord$TS), ]
# 
# ETr_filt_FILE <- ETr_smth_IP
# subs.d <- as.matrix(t(ETr_filt_ord[,-1]))
# 
# ETr_filt_FILE[9:nrow(ETr_filt_FILE), 7:ncol(ETr_filt_FILE)] <- subs.d
# 
# write.csv(ETr_filt_FILE, paste0(opPATH, "OP-18","_ETr_Thresh_Filt.csv"))
# 
# 
# # 5. Thresholded/Filtered ETr interpolation
# ETr_filt_imputed_FILE <- ETr_filt_FILE
# subs.d.imp <- subs.d
# na.sum <- c()
# 
# for(i in 1:nrow(subs.d.imp))
# {
#   subs.d.imp[i, ] <- na.aggregate.default(subs.d.imp[i, ])
#   na.sum[i] <- sum(is.na(subs.d.imp[i, ]))
# }
# print(paste0("#Sectors still with NA: ", (length(na.sum) - length(!na.sum == 0))))
# 
# ETr_filt_imputed_FILE[9:nrow(ETr_filt_imputed_FILE), 7:ncol(ETr_filt_imputed_FILE)] <- subs.d.imp
# 
# write.csv(ETr_filt_imputed_FILE, paste0(opPATH, "OP-19","_ETr_Thresh_Impt.csv"))

###########################################################

##### Start Tr extraction to get raw Tr profiles #####
# pe.df$TS <- dmy_hm(pe.df$timestamp) # Check format and set TS
# 
# pe.df$date <- date(pe.df$TS)
# 
# sel.secs <- unique(na.omit(ETr_filt_imputed_FILE$old_unit))
# 
# pe.df.ETr <- pe.df[pe.df$date %in% unq.dts & pe.df$Sector %in% sel.secs, ]
# 
# # Select only the req. Treatment/Mulch type, here Mulch2
# pe.df.ETr <- pe.df.ETr[pe.df.ETr$Treatment=="Mulch2", ]
# sel.secs <- unique(pe.df.ETr$Sector) # update after selection
# 
# # Select the following columns: "Sector", "Treatment", "Genotype", "LeafArea3D", "TS", "date"
# pe.df.ETr <- pe.df.ETr[ ,c(1, 3, 5, 12, 17, 18)]
# 
# colnames(pe.df.ETr) <- c("old_unit", "Treatment", "Genotype", "LeafArea3D", "TS", "date")
# pe.df.ETr$Genotype <- factor(pe.df.ETr$Genotype)
# pe.df.ETr$Treatment <- factor(pe.df.ETr$Treatment)
# 
# pe.ETr.grpDT <- pe.df.ETr %>% group_by(old_unit, date, Treatment, Genotype) %>% 
#                 dplyr::summarise(Max = max(LeafArea3D, na.rm=TRUE))
# 
# names(pe.ETr.grpDT)[5] <- "LeafArea3D"
# 
# # make a matrix of LA3D of dim = ETr_core data, then calculate TR #
# ## LAI=((((3DLA/100) - 29.9)/0.36)*(1/0.26))/10000 ##
# ## T = (1-exp(-0.463*LAI))*ETr ##
# 
# LAI.mat <- matrix(NA, nrow = length(sel.secs), ncol = ncol(subs.d))
# rownames(LAI.mat) <- sel.secs
# 
# # sort LAI.mat rownames as per the ETR_smooth file 
# LAI.mat <- LAI.mat[order(match(rownames(LAI.mat), 
#                   ETr_filt_imputed_FILE[9:nrow(ETr_filt_imputed_FILE), "old_unit"])),]
# 
# LAI.all.dates <- LAI.mat; unq.dts.copy <- unq.dts
# 
# LAI.all.dates <- LAI.all.dates[ , c(1:length(unq.dts))]
# colnames(LAI.all.dates) <- c(as.character(unq.dts.copy))
# 
# write.csv(ETr_filt_imputed_FILE, "OP-20","SG_Tr_IP.csv")
# 
# # Calculate raw Transpiration, Tr
# 
# # Select the desired "Treatment" type e.g. "Mulch2"
# M2rows <- which(ETr_filt_imputed_FILE$Treatment=="Mulch2") 
# Tr_IP_Mulch2 <- ETr_filt_imputed_FILE[c(1:8, M2rows), ]
# 
# Tr_OP <- calculateTr(x = Tr_IP_Mulch2, y = pe.df.ETr, 
#                      z = LAI.mat, d = unq.dts) #If >1 Treatment type,
#                      # else x = ETr_filt_imputed_FILE
# 
# raw.trans.mat <- Tr_OP$Trans.mat
# 
# LA3D.all.dates <- Tr_OP$LA3D_TS
# 
# LAI.mat <- Tr_OP$LAI.mat
# 
# 
# raw.trans <- Tr_IP_Mulch2 #else raw.trans = ETr_filt_imputed_FILE
# 
# raw.trans[9:nrow(raw.trans), 7:ncol(raw.trans)] <- raw.trans.mat
# 
# write.csv(raw.trans, paste0(opPATH, "OP-21","_raw_Tr.csv"))
# 
# write.csv(LA3D.all.dates, paste0(opPATH, "OP-22","_3DLA_TS.csv"))
# 
# write.csv(LAI.mat, paste0(opPATH, "OP-23","_LAI_TS.csv"))


##### Feature Extraction of RAW Transpiration Data #####
# featuresRES <- getFeatures(x = raw.trans)
featuresRES <- getFeatures(x = wthr.ETref.ETobs)

allFeatures <- featuresRES$allFeatures
unq.dts <- unique(substr(colnames(wthr.ETref.ETobs)[-c(1:6)], 1, 10))

# create H2 dataframe to store H2 est. of each feature for each day
F.He <- as.data.frame(matrix(NA, nrow = length(unq.dts), ncol = 15)) # Date-ROW, feature-COL
colnames(F.He) <- c("maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
                    "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
                    "auc.prop.07-19", "auc.night", "cos.sim.index")
rownames(F.He) <- unq.dts

featureHeRES <- getFeatureHe(x = allFeatures, y = wthr.ETref.ETobs, d = unq.dts, p =opPATH.raw)
write.csv(featureHeRES, paste0(opPATH, "rawETr_featureH2.csv"))

### save all features as feature Time Series ###
### Each feature set: dim(length(unq.dts) x (nrow(raw.trans)-8)) ###

## Prepare data for 'each feature'
maxET <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
slope.maxET.6 <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
slope.07maxET <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
slope.00.07 <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
slope.19.2345 <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
curvmaxET <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
total.auc <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
auc.10.15 <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
sd.10.15 <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
auc.prop.10.15 <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
auc.07.19 <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
sd.07.19 <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
auc.prop.07.19 <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
auc.night <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))
cos.sim.index <- as.data.frame(matrix(nr = (nrow(wthr.ETref.ETobs)-8), nc = length(unq.dts)))

for (j in 1:(nrow(wthr.ETref.ETobs)-8)){
  
  # "maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
  # "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
  # "auc.prop.07-19", "auc.night", "cos.sim.index"
  
  for(d in 1:length(unq.dts))
  {maxET[j, d] <- data.frame(allFeatures[[j]][d, 1])
  slope.maxET.6[j, d] <- data.frame(allFeatures[[j]][d, 2])
  slope.07maxET[j, d] <- data.frame(allFeatures[[j]][d, 3])
  slope.00.07 [j, d] <- data.frame(allFeatures[[j]][d, 4])
  slope.19.2345[j, d] <- data.frame(allFeatures[[j]][d, 5])
  
  curvmaxET[j, d] <- data.frame(allFeatures[[j]][d, 6])
  total.auc[j, d] <- data.frame(allFeatures[[j]][d, 7])
  auc.10.15[j, d] <- data.frame(allFeatures[[j]][d, 8])
  sd.10.15 [j, d] <- data.frame(allFeatures[[j]][d, 9])
  auc.prop.10.15[j, d] <- data.frame(allFeatures[[j]][d, 10])
  
  auc.07.19[j, d] <- data.frame(allFeatures[[j]][d, 11])
  sd.07.19[j, d] <- data.frame(allFeatures[[j]][d, 12])
  auc.prop.07.19[j, d] <- data.frame(allFeatures[[j]][d, 13])
  auc.night [j, d] <- data.frame(allFeatures[[j]][d, 14])
  cos.sim.index[j, d] <- data.frame(allFeatures[[j]][d, 15])
  }
  
} 

names(maxET)=names(slope.maxET.6)=names(slope.07maxET)=names(slope.00.07)=names(slope.19.2345)<-unq.dts
names(curvmaxET)=names(total.auc)=names(auc.10.15)=names(sd.10.15)=names(auc.prop.10.15)<-unq.dts
names(auc.07.19)=names(sd.07.19)=names(auc.prop.07.19)=names(auc.night)=names(cos.sim.index)<-unq.dts


write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], maxET)), 
          paste0(opPATH.raw, "maxET.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], slope.maxET.6)), 
          paste0(opPATH.raw, "slope.maxET.6.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], slope.00.07)), 
          paste0(opPATH.raw, "slope.00.07.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], slope.07maxET)), 
          paste0(opPATH.raw, "slope.07maxET.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], slope.19.2345)), 
          paste0(opPATH.raw, "slope.19.2345.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], curvmaxET)), 
          paste0(opPATH.raw, "curvmaxET.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], total.auc)), 
          paste0(opPATH.raw, "total.auc.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], auc.10.15)), 
          paste0(opPATH.raw, "auc.10.15.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], sd.10.15)), 
          paste0(opPATH.raw, "sd.10.15.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], auc.prop.10.15)), 
          paste0(opPATH.raw, "auc.prop.10.15.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], auc.07.19)), 
          paste0(opPATH.raw, "auc.07.19.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], sd.07.19)), 
          paste0(opPATH.raw, "sd.07.19.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], auc.prop.07.19)), 
          paste0(opPATH.raw, "auc.prop.07.19.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], auc.night)), 
          paste0(opPATH.raw, "auc.night.csv"))

write.csv(as.data.frame(cbind(wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs), 1:6], cos.sim.index)), 
          paste0(opPATH.raw, "cos.sim.index.csv"))



###### Stage III: Generate smooth ETr and obtain smooth Tr features ######
# ETr smoothing
ETr_smthIP <- wthr.ETref.ETobs

subs.d.imp <- ETr_smthIP[9:nrow(ETr_smthIP),
                          7:ncol(ETr_smthIP)]

smoothETrMAT <- smoothETr(x = subs.d.imp)


ETr_smoothFILE <- ETr_smthIP
ETr_smoothFILE[9:nrow(ETr_smoothFILE), 
               7:ncol(ETr_smoothFILE)] <- smoothETrMAT

write.csv(ETr_smoothFILE, paste0(opPATH, "OP-9","_ETr_smth.csv"))


# Calculate Tr from smooth ETr
# Tr_OP <- calculateTr(x = ETr_smoothFILE, y = pe.df.ETr, z = LAI.mat, d = unq.dts)
# Tr_OP <- smoothETr(x = raw.trans.mat)
# 
# smth.trans <- raw.trans
# 
# smth.trans[9:nrow(smth.trans), 7:ncol(smth.trans)] <- Tr_OP/1000
# 
# write.csv(smth.trans, paste0(opPATH, "OP-25","_smth_Tr.csv"))


### Feature Extraction of SMOOTH Transpiration Data ###
featuresRES <- getFeatures(x = ETr_smoothFILE)

allFeatures <- featuresRES$allFeatures

# create H2 dataframe to store H2 est. of each feature for each day
F.He <- as.data.frame(matrix(NA, nrow = length(unq.dts), ncol = 15)) # Date-ROW, feature-COL
colnames(F.He) <- c("maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
                    "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
                    "auc.prop.07-19", "auc.night", "cos.sim.index")
rownames(F.He) <- unq.dts

featureHeRES <- getFeatureHe(x = allFeatures, y = ETr_smoothFILE, d = unq.dts, p = opPATH.smth)
write.csv(featureHeRES, paste0(opPATH, "smthETr_featureH2.csv"))


### save all features as feature Time Series ###
### Each feature set: dim(length(unq.dts) x (nrow(raw.trans)-8)) ###

## Prepare data for 'each feature'
maxET <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
slope.maxET.6 <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
slope.07maxET <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
slope.00.07 <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
slope.19.2345 <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
curvmaxET <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
total.auc <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
auc.10.15 <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
sd.10.15 <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
auc.prop.10.15 <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
auc.07.19 <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
sd.07.19 <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
auc.prop.07.19 <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
auc.night <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))
cos.sim.index <- as.data.frame(matrix(nr = (nrow(ETr_smoothFILE)-8), nc = length(unq.dts)))

for (j in 1:(nrow(ETr_smoothFILE)-8)){
  
  # "maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
  # "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
  # "auc.prop.07-19", "auc.night", "cos.sim.index"
  
  for(d in 1:length(unq.dts))
  {maxET[j, d] <- data.frame(allFeatures[[j]][d, 1])
  slope.maxET.6[j, d] <- data.frame(allFeatures[[j]][d, 2])
  slope.07maxET[j, d] <- data.frame(allFeatures[[j]][d, 3])
  slope.00.07 [j, d] <- data.frame(allFeatures[[j]][d, 4])
  slope.19.2345[j, d] <- data.frame(allFeatures[[j]][d, 5])
  
  curvmaxET[j, d] <- data.frame(allFeatures[[j]][d, 6])
  total.auc[j, d] <- data.frame(allFeatures[[j]][d, 7])
  auc.10.15[j, d] <- data.frame(allFeatures[[j]][d, 8])
  sd.10.15 [j, d] <- data.frame(allFeatures[[j]][d, 9])
  auc.prop.10.15[j, d] <- data.frame(allFeatures[[j]][d, 10])
  
  auc.07.19[j, d] <- data.frame(allFeatures[[j]][d, 11])
  sd.07.19[j, d] <- data.frame(allFeatures[[j]][d, 12])
  auc.prop.07.19[j, d] <- data.frame(allFeatures[[j]][d, 13])
  auc.night [j, d] <- data.frame(allFeatures[[j]][d, 14])
  cos.sim.index[j, d] <- data.frame(allFeatures[[j]][d, 15])
  }
  
} 

names(maxET)=names(slope.maxET.6)=names(slope.07maxET)=names(slope.00.07)=names(slope.19.2345)<-unq.dts
names(curvmaxET)=names(total.auc)=names(auc.10.15)=names(sd.10.15)=names(auc.prop.10.15)<-unq.dts
names(auc.07.19)=names(sd.07.19)=names(auc.prop.07.19)=names(auc.night)=names(cos.sim.index)<-unq.dts


write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], maxET)), 
          paste0(opPATH.smth, "maxET.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], slope.maxET.6)), 
          paste0(opPATH.smth, "slope.maxET.6.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], slope.00.07)), 
          paste0(opPATH.smth, "slope.00.07.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], slope.07maxET)), 
          paste0(opPATH.smth, "slope.07maxET.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], slope.19.2345)), 
          paste0(opPATH.smth, "slope.19.2345.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], curvmaxET)), 
          paste0(opPATH.smth, "curvmaxET.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], total.auc)), 
          paste0(opPATH.smth, "total.auc.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], auc.10.15)), 
          paste0(opPATH.smth, "auc.10.15.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], sd.10.15)), 
          paste0(opPATH.smth, "sd.10.15.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], auc.prop.10.15)), 
          paste0(opPATH.smth, "auc.prop.10.15.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], auc.07.19)), 
          paste0(opPATH.smth, "auc.07.19.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], sd.07.19)), 
          paste0(opPATH.smth, "sd.07.19.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], auc.prop.07.19)), 
          paste0(opPATH.smth, "auc.prop.07.19.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], auc.night)), 
          paste0(opPATH.smth, "auc.night.csv"))

write.csv(as.data.frame(cbind(ETr_smoothFILE[9:nrow(ETr_smoothFILE), 1:6], cos.sim.index)), 
          paste0(opPATH.smth, "cos.sim.index.csv"))

##### Smooth Transpiration Rate (TR) Calculation #####

ets = ETr_smoothFILE

laidf = read.csv("./data/PERIODE2-01-10-10-10/LA_PERIODE2.csv", sep = ";")
#laidf <- na.approx(laidf[,9:ncol(laidf)])
# Get unique sector names
ets.sec <- unique(ets$old_unit)
lai_sec <- laidf$old_unit

# Convert column names to Date format
colnames(laidf)[9:ncol(laidf)] <- sub("X","",colnames(laidf)[9:ncol(laidf)])

lai_dts <- colnames(laidf)[9:ncol(laidf)]
lai_dtsseq <- dmy_hm(lai_dts)

colnames(laidf)[9:ncol(laidf)]<-as.character(lai_dtsseq) 

# Extract ETr part
etsMat <- ets[9:nrow(ets), c(1,7:ncol(ets))]
trMat <- etsMat # TR placeholder matrix

# TR calculation
for(i in 1:nrow(etsMat)) {
  print(i)
  etsTmp <- etsMat[i,-1]
  etsTmpMat <- matrix(etsTmp, nrow = length(etsTmp)/length(unq.dts))
  
  for(j in 1:length(unq.dts)){
    print(j)
    d_j <- unq.dts[j]
    sec_j <- etsMat[i,1]
    lai_j <- laidf[laidf$unit==sec_j, d_j]
    
    if(is.na(lai_j)==TRUE){
      lai_j == 0.0
      etsTmpMat[,j] <- unlist(etsTmpMat[,j])/lai_j
    }else{
      etsTmpMat[,j] <- unlist(etsTmpMat[,j])/lai_j
    }
  }
  
  trMat[i,2:ncol(trMat)]<-c(etsTmpMat)
}

TR_smoothFILE <- ETr_smoothFILE
TR_smoothFILE[9:nrow(TR_smoothFILE), 
              7:ncol(TR_smoothFILE)] <- trMat[ , 2:ncol(trMat)]

write.csv(TR_smoothFILE, paste0(opPATH, "OP-10_TR_smth.csv"))

### Feature Extraction of SMOOTH TR Data ###
featuresRES <- getFeatures(x = TR_smoothFILE)

allFeatures <- featuresRES$allFeatures

# create H2 dataframe to store H2 est. of each feature for each day
F.He <- as.data.frame(matrix(NA, nrow = length(unq.dts), ncol = 15)) # Date-ROW, feature-COL
colnames(F.He) <- c("maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
                    "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
                    "auc.prop.07-19", "auc.night", "cos.sim.index")
rownames(F.He) <- unq.dts

featureHeRES <- getFeatureHe(x = allFeatures, y = ETr_smoothFILE, 
                             d = unq.dts, p = opPATH.smthTR)
write.csv(featureHeRES, paste0(opPATH, "smthTR_featureH2.csv"))


### save all features as feature Time Series ###
### Each feature set: dim(length(unq.dts) x (nrow(raw.trans)-8)) ###

## Prepare data for 'each feature'
maxET <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
slope.maxET.6 <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
slope.07maxET <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
slope.00.07 <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
slope.19.2345 <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
curvmaxET <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
total.auc <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
auc.10.15 <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
sd.10.15 <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
auc.prop.10.15 <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
auc.07.19 <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
sd.07.19 <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
auc.prop.07.19 <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
auc.night <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))
cos.sim.index <- as.data.frame(matrix(nr = (nrow(TR_smoothFILE)-8), nc = length(unq.dts)))

for (j in 1:(nrow(TR_smoothFILE)-8)){
  
  # "maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
  # "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
  # "auc.prop.07-19", "auc.night", "cos.sim.index"
  
  for(d in 1:length(unq.dts))
  {maxET[j, d] <- data.frame(allFeatures[[j]][d, 1])
  slope.maxET.6[j, d] <- data.frame(allFeatures[[j]][d, 2])
  slope.07maxET[j, d] <- data.frame(allFeatures[[j]][d, 3])
  slope.00.07 [j, d] <- data.frame(allFeatures[[j]][d, 4])
  slope.19.2345[j, d] <- data.frame(allFeatures[[j]][d, 5])
  
  curvmaxET[j, d] <- data.frame(allFeatures[[j]][d, 6])
  total.auc[j, d] <- data.frame(allFeatures[[j]][d, 7])
  auc.10.15[j, d] <- data.frame(allFeatures[[j]][d, 8])
  sd.10.15 [j, d] <- data.frame(allFeatures[[j]][d, 9])
  auc.prop.10.15[j, d] <- data.frame(allFeatures[[j]][d, 10])
  
  auc.07.19[j, d] <- data.frame(allFeatures[[j]][d, 11])
  sd.07.19[j, d] <- data.frame(allFeatures[[j]][d, 12])
  auc.prop.07.19[j, d] <- data.frame(allFeatures[[j]][d, 13])
  auc.night [j, d] <- data.frame(allFeatures[[j]][d, 14])
  cos.sim.index[j, d] <- data.frame(allFeatures[[j]][d, 15])
  }
  
} 

names(maxET)=names(slope.maxET.6)=names(slope.07maxET)=names(slope.00.07)=names(slope.19.2345)<-unq.dts
names(curvmaxET)=names(total.auc)=names(auc.10.15)=names(sd.10.15)=names(auc.prop.10.15)<-unq.dts
names(auc.07.19)=names(sd.07.19)=names(auc.prop.07.19)=names(auc.night)=names(cos.sim.index)<-unq.dts


write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], maxET)), 
          paste0(opPATH.smthTR, "maxET.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], slope.maxET.6)), 
          paste0(opPATH.smthTR, "slope.maxET.6.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], slope.00.07)), 
          paste0(opPATH.smthTR, "slope.00.07.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], slope.07maxET)), 
          paste0(opPATH.smthTR, "slope.07maxET.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], slope.19.2345)), 
          paste0(opPATH.smthTR, "slope.19.2345.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], curvmaxET)), 
          paste0(opPATH.smthTR, "curvmaxET.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], total.auc)), 
          paste0(opPATH.smthTR, "total.auc.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], auc.10.15)), 
          paste0(opPATH.smthTR, "auc.10.15.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], sd.10.15)), 
          paste0(opPATH.smthTR, "sd.10.15.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], auc.prop.10.15)), 
          paste0(opPATH.smthTR, "auc.prop.10.15.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], auc.07.19)), 
          paste0(opPATH.smthTR, "auc.07.19.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], sd.07.19)), 
          paste0(opPATH.smthTR, "sd.07.19.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], auc.prop.07.19)), 
          paste0(opPATH.smthTR, "auc.prop.07.19.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], auc.night)), 
          paste0(opPATH.smthTR, "auc.night.csv"))

write.csv(as.data.frame(cbind(TR_smoothFILE[9:nrow(TR_smoothFILE), 1:6], cos.sim.index)), 
          paste0(opPATH.smthTR, "cos.sim.index.csv"))

end.time <- Sys.time()

print(paste0("Complete processing executed in: ", round((end.time-st.time), 2), "minutes"))

