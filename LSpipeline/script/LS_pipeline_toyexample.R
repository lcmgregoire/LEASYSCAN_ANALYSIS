############################################
# Overall LeasyScan pipeline - toy example #
############################################

library("devtools")

#install_github("ICRISAT-GEMS/LoadCellDataProcessing", force = TRUE)
#install_github("ICRISAT-GEMS/platformDataAnalysis", force = TRUE)

library("dplyr")
library("LoadCellDataProcessing")
library("statgenHTP")
library("platformDataAnalysis")
library("ggplot2")
library("lubridate")

# define reference WD
ref_wd <- "C:/Users/2021lg003/Documents/Sorghum/Experiments/EXP_INDE_2023/ANALYSIS_ALL_DATASET/LSpipeline/data"

# data preparation ----

# PE data
setwd('C:/Users/2021lg003/Documents/Sorghum/Experiments/EXP_INDE_2023/DATASET/FULL_DATASET')
pe_data <- read.csv(file = 'Exp60 Sorghum Ref set IRD Trial Sep 2023-477 mm_20231027_planteye.csv', sep = ";")
pe_data<-pe_data[,-(1:4)] 


# LC data 
setwd('C:/Users/2021lg003/Documents/Sorghum/Experiments/EXP_INDE_2023/DATASET/FULL_DATASET')
lc_data <- read.csv(file = 'Exp60 Sorghum Ref set IRD Trial Sep 2023-477 mm_20231027_DroughtSpotter.csv', sep = ";")
#!! different from template !
##format should be : unit, sensor, variable, timestamp, value (all chr except value) 
## !! careful with timestamp : yyyy-mm-dd hh:mm in character 
lc_data<-lc_data[,-(1:4)]
lc_data<-cbind(lc_data$unit, lc_data$sensor, "Weight g", lc_data$timestamp, lc_data$Weight.g)
lc_data<-as.data.frame(lc_data)
colnames(lc_data)<-c("unit", "sensor", "variable", "timestamp", "value")
lc_data$value<-as.numeric(lc_data$value)

# exp design 
setwd('C:/Users/2021lg003/Documents/Sorghum/Experiments/EXP_INDE_2023/DATASET/FULL_DATASET')

d_design <- read.csv(file = "Exp60_design_info.csv", sep = ";")
d_design <- d_design %>% select(row.psx, col.psx, rep, block, geno_id)
d_exp<-d_design

# select a small subset
#d_des_s <- d_design %>% filter(row.psx %in% 3:6, col.psx %in% 1:12)
#d_des_s <- d_design %>% filter(row.psx %in% 1:3, col.psx %in% 1:12) ##work on the row 1-3 , the whole column 
#d_exp <- d_des_s

# save data
setwd(ref_wd)
save(d_exp, file = "exp_des.RData")
## to see how it looks like : 
load(file = "exp_des.RData")

# convert the row and columns in reference system: col.psx -> row
# row.psx -> col + 2 (the indice depend from where we start here A1 and A2 are skipped)

d_exp$rowNum <- d_exp$col.psx
d_exp$colNum <- d_exp$row.psx + 2
d_exp$plotId <- paste0(paste0("c", d_exp$colNum), paste0("r", d_exp$rowNum))

# add the sector id to the experimental design
setwd("C:/Users/2021lg003/Documents/Sorghum/Experiments/EXP_INDE_2023/ANALYSIS_ALL_DATASET/LSpipeline/data")
load(file = "unit_row_col_block.RData")
#d_unit=read.csv("sector.csv", sep = ";")
d_unit <- d_unit_row_col_block[, c(2:4, 8)]
d_exp <- left_join(x = d_exp, y =  d_unit, by = "plotId")
#d_exp$new_unit= d_exp$new_unit.x ### if duplication of new.unit column

# subset pe, lc, weather, and sensor data
pe_data_s <- pe_data[pe_data$unit %in% d_exp$new_unit, ]

# check genotypes correspond to the experimental design
substr(unique(pe_data_s$genotype), 11, nchar(unique(pe_data_s$genotype)))

pe_data <- pe_data_s

lc_data_s <- lc_data[lc_data$unit %in% d_exp$new_unit, ]

unique(lc_data_s$unit)

lc_data <- lc_data_s

# weather data
setwd('C:/Users/2021lg003/Documents/Sorghum/Experiments/EXP_INDE_2023/DATASET/FULL_DATASET')
wth_data <- read.csv(file = 'Exp60 Sorghum Ref set IRD Trial Sep 2023-477 mm_20231027_Climate Datalogger.csv')

# sensor unit map
setwd('C:/Users/2021lg003/Documents/Sorghum/Experiments/EXP_INDE_2023/DATASET/FULL_DATASET')
sensor_data <- read.csv(file = 'Exp60 Sorghum Ref set IRD Trial Sep 2023-477 mm_20231027_sensorUnitMap.csv')

# select a subset of the sensor data through unit present in the design
sensor_data <- sensor_data[sensor_data$unit %in% d_exp$new_unit, ]

# subset the weather data given the unique sensor
wth_data <- wth_data[wth_data$sensor %in% unique(sensor_data$sensor), ]

# save data
setwd(ref_wd)
save(pe_data, file = "pe_data.RData")
save(lc_data, file = "lc_data.RData")
save(sensor_data, file = "sensor_data.RData")
save(wth_data, file = "weather_data.RData")
#######################################################" COMPLETED 
# Initial filtering through experimental design ----

# load data
setwd(ref_wd)
load(file = "exp_des.RData")
load(file = "pe_data.RData")
load(file = "lc_data.RData")
load(file = "sensor_data.RData")
load(file = "weather_data.RData")

# starting point: experimental design: list of genotype, positions (row, col), ...
head(d_exp)

# need to add unit (sector) identifiers if not present in the design sheet

# convert the row and columns in reference system: col.psx -> row
# row.psx -> col + 2 (the indice depend from where we start here A1 and A2 are skipped)

d_exp$rowNum <- d_exp$col.psx
d_exp$colNum <- d_exp$row.psx + 2
d_exp$plotId <- paste0(paste0("c", d_exp$colNum), paste0("r", d_exp$rowNum))

# reference data with conversion sector id given row and column position
load(file = "unit_row_col_block.RData")

# GET row col from sector (unit) identifiers

# d_exp_test <- data.frame(unit = d_unit_row_col_block$new_unit[1:10],
#                          geno = paste0("G", 1:10))
# # ajouter row and col info
# d_exp_test$rowNum <- LS_new_unit_row_lk[d_exp_test$unit]
# d_exp_test$colNum <- LS_new_unit_col_lk[d_exp_test$unit]

d_unit <- d_unit_row_col_block[, c(2:4, 8)]
d_exp <- left_join(x = d_exp, y =  d_unit, by = "plotId")

# filter pe_data using exp_des unit information
pe_data <- pe_data[pe_data$unit %in% d_exp$new_unit, ]

# format the genotype data (match further with other data)
genotype <- strsplit(x = pe_data$genotype, split = '_')
geno_n_pieces <- unlist(lapply(X = genotype, FUN = length))
table(geno_n_pieces)

geno_vec <- unlist(lapply(X = genotype, FUN = `[[`, 3))
unique(geno_vec)
table(geno_vec)

# geno_vec[geno_vec == ""] <- NA
pe_data$genotype <- geno_vec

# review the genotype name check using d_exp as reference and link through unit.

# filter lc_data using exp_des unit information
lc_data <- lc_data[lc_data$unit %in% d_exp$new_unit, ] ## CAREFUL WITH THE 0 in LC unit. Ex: 131:01:02 or 131:1:2 

# add a genotype column using exp design sector genotype look_up
geno_lk <- d_exp$geno_id
names(geno_lk) <- d_exp$new_unit
lc_data$genotype <- lc_data$g_alias <- geno_lk[lc_data$unit]

# select a subset of the sensor data through unit present in the design
sensor_data <- sensor_data[sensor_data$unit %in% d_exp$new_unit, ]

# subset the weather data given the unique sensor
wth_data <- wth_data[wth_data$sensor %in% unique(sensor_data$sensor), ]

# save the data
setwd(ref_wd)
save(d_exp, file = "exp_des_filtered.RData")
save(pe_data, file = "pe_data_filtered.RData")
save(lc_data, file = "lc_data_filtered.RData")
save(sensor_data, file = "sensor_data_filtered.RData")
save(wth_data, file = "weather_data_filtered.RData")
####################################################### COMPLETED
############################ LC pipeline #######################################
# LC pipeline: Data processing ----

setwd(ref_wd)
load(file = "pe_data_filtered.RData")
load(file = "lc_data_filtered.RData")
load(file = "sensor_data_filtered.RData")
load(file = "weather_data_filtered.RData")

# PE_data
pe_data <- pe_data %>% select(unit, genotype, timestamp, Leaf.area..mmÂ.. ) ## or Leaf.area..mm..
pe_data$replicate <- NA
colnames(pe_data) <- c('sector', 'genotype', 'timestamp', 'leaf_area', 'replicate')
col_names <- c("sector", "genotype", "replicate", "timestamp", "leaf_area")
pe_data <- pe_data[, col_names]


# LC data 
lc_data <- lc_data[lc_data$variable == "Weight g", ] ##We have only weight here ! nothing happens
lc_data$treatment <- NA
lc_data <- lc_data %>% rename(`Mass..g` = value, sector = unit)
lc_col_name <- c("sector", "genotype", "g_alias", "treatment", "timestamp", "Mass..g")
lc_data <- lc_data[, lc_col_name]


# sensor data #required format: sector, sensor
sensor_data <- sensor_data[, 1:2]
colnames(sensor_data) <- c("sector", "sensor")


# weather data

# modify few climatic variable name to be exactly as in the example
wth_data$variable[wth_data$variable == "Relative Humidity (%)"] <- "Relative humidity (%)"
wth_data$variable[wth_data$variable == "Solar Radiation (W/(s*m2))"] <- "Solar radiation (W/(s*m²))"
wth_data$variable[wth_data$variable == "Wind Direction (°)"] <- "Wind direction (°)"


unique(wth_data$variable)

# save the data
setwd(ref_wd)
save(pe_data, file = "pe_data_processed.RData")
save(lc_data, file = "lc_data_processed.RData")
save(sensor_data, file = "sensor_data_processed.RData")
save(wth_data, file = "weather_data_processed.RData")
####################################COMPLETED §################
# LC pipeline: TR feature data extraction (TR_data_proc) ----

# load data
setwd(ref_wd)
load(file = "pe_data_processed.RData")
load(file = "lc_data_processed.RData")
load(file = "sensor_data_processed.RData")
load(file = "weather_data_processed.RData")

#Error in check_timestamp_format(d$timestamp) : The date in the vector do not all match the %Y-%m-%d %H:%M:%S format.
#CAREFUL WHEN YOU LOAD THE DATA ! change in excel for a ymd_hms format for all files. 


t1 <- Sys.time()

TR_res <- TR_data_proc(lc_data = lc_data, pe_data = pe_data,
                       wth_data = wth_data, sensor_data = sensor_data,
                       lastDate = NULL, skew_test = FALSE,
                       LAI_correction = FALSE)

t2 <- Sys.time()

t_diff <- t2 - t1
print(t_diff)

# save the results
save(TR_res, file = "TR_res.RData")

# LC pipeline: check the results shape, possible selection of days ----

# load data
setwd(ref_wd)
load(file = "TR_res.RData")

p <- plot_whole_TR_time_series(results = TR_res, sector_sel = c(2, 4, 8))
p

# loop over all sectors to make a pdf with all plots
n_sect <- nrow(TR_res$TR_smth$TR_smth)

pdf(file = "./plot/TR_time_series_full.pdf", width = 12, height = 8)

for(i in 1:n_sect){

  p <- tryCatch(plot_whole_TR_time_series(results = TR_res, sector_sel = i,
                                          main = paste("sector", i)),
                error = function(e) NULL)

  if(!is.null(p)){ print(p) }

}

dev.off()

# leaf area index coverage
# p <- plot_LAI_coverage(results = TR_res)
# p

# check Tr feature daily TS
trait_id <- "TRmax"
title_i <- paste(trait_id, "- KK subset")

p <- plot_TR_time_series(results = TR_res, trait = trait_id, n_sector = NULL,
                         color = FALSE, main = title_i)
p


# select days
start_day = 1
end_day = 8
TR_feature <- TR_res$TR_smth$Max_TR_smth
ts_date <- colnames(TR_feature)[6:ncol(TR_feature)]
LC_day_sel <- ts_date[start_day:end_day]

# the days can be used later to subset
TR_res <- TR_res_subset(TR_res = TR_res, start_day = start_day, end_day = end_day)

# save the results
save(TR_res, file = "TR_res_sub.RData")

# LC pipeline: TR results processing TP object ----

setwd(ref_wd)
load(file = "exp_des_filtered.RData")
load(file = "TR_res_sub.RData")

d_TP <- TR_res_to_pre_TP(TR_res)

# extract the TR ~ VPD regression coefficient
TR_VPD_reg_coeff <- TR_VPD_reg(TR_res = TR_res, do_plot = FALSE)
TR_VPD_lk <- TR_VPD_reg_coeff$Tr_VPD_slope
names(TR_VPD_lk) <- TR_VPD_reg_coeff$unit
d_TP$TR_VPD <- TR_VPD_lk[d_TP$unit]

# add time number
d_TP$timestamp <- strptime(d_TP$timestamp, "%Y-%m-%d")
d_TP <- add_timeNumber(d_TP)
d_TP <- d_TP %>% arrange(unit, timestamp)

# add extra columns from the experimental design
d_TP <- add_exp_des_col(data = d_TP, d_exp_des = d_exp,
                        data_unit = "unit",
                        d_exp_unit = "new_unit",
                        col_add = c("rowNum", "colNum", "block", "cross"))

# PlotId is the combination of row and col information
d_TP$plotId <- paste0(paste0("c", d_TP$colNum), paste0("r", d_TP$rowNum))

# arrange the columns in a certain order
d_TP_meta <- d_TP %>% select(timeNumber, timestamp, block, rowNum, colNum, plotId, cross,
                             genotype) %>% rename(timePoint = timestamp)

d_TP_trait <- d_TP %>% select(contains("_smth"), TR_VPD)

# modify the names of the traits
tr_nm <- c("max_TR", "slope_6pt_bfr_maxTR", "slope_07_maxTR",
           "slope_00_07", "slope_19h_23h45", "curve_maxTR",
           "total_auc", "auc_10h_15h", "sd_10h_15h",
           "prop_auc_10h_15h", "auc_7h_19h", "sd_7h_19h",
           "prop_auc_7h_19h", "auc_night", "TR_VPD_slope")

colnames(d_TP_trait) <- tr_nm

data_LC <- data.frame(d_TP_meta, d_TP_trait)

TP_LC <- createTimePoints(dat = data_LC,
                          experimentName = "Exp_51_Kenin_Keni_LC",
                          genotype = "genotype",
                          timePoint = "timePoint",
                          plotId = "plotId",
                          rowNum = "rowNum", colNum = "colNum")

# LC pipeline: diagnostic plot ----

plot(TP_LC, traits = "max_TR",
     plotType = "layout",
     timePoints = 1)

# LC pipeline: TR results spatial adjustment ----

LC_adjusted <- spatial_adjustment(TP = TP_LC,
                                  traits = c("max_TR", "total_auc", "TR_VPD_slope"),
                                  timePoints = NULL,
                                  extraFixedFactors = "block",
                                  geno.decomp = NULL,
                                  what = "fixed",
                                  useCheck = FALSE,
                                  useRepId = FALSE,
                                  engine = "SpATS",
                                  spatial = FALSE,
                                  quiet = TRUE)

plot_LC_adj <- LC_adjusted$plot_res
geno_LC_adj <- LC_adjusted$geno_res
comp_mon <- LC_adjusted$comp_monitor
h2_LC <- LC_adjusted$h2_res

# check adjusted TS

plot_trend(data = plot_LC_adj, trait = "max_TR_corr")
plot_trend(data = plot_LC_adj, trait = "total_auc_corr")

plot_trend(data = geno_LC_adj, trait = "max_TR_pred", genotype = TRUE)
plot_trend(data = geno_LC_adj, trait = "total_auc_pred", genotype = TRUE)

# save the data
setwd(ref_wd)
save(plot_LC_adj, file = "LC_plot_adjusted.RData")
save(geno_LC_adj, file = "LC_geno_adjusted.RData")
save(h2_LC, file = "h2_LC.RData")


############################ PE pipeline #######################################
# PE pipeline: Data processing ----

setwd(ref_wd)
load(file = "pe_data_filtered.RData")
load(file = "exp_des_filtered.RData")

colnames(pe_data) <- mdf_raw_pe_colnames(colnames = colnames(pe_data))

ref_trait_nm <- c("Digital_biomass", "Height", "Leaf_angle", "Leaf_area",
                  "Leaf_area_index", "Leaf_area_projected", "Leaf_inclination",
                  "Light_penetration_depth")

pe_data <- pe_data %>% select(unit, genotype, g_alias, treatment, timestamp,
                              Digital_biomass, Height, Leaf_angle, Leaf_area, Leaf_area_index,
                              Leaf_area_projected, Leaf_inclination, Light_penetration_depth)

colnames(pe_data)[6:13] <- ref_trait_nm

# add time number
pe_data$timestamp <- as.POSIXlt(pe_data$timestamp, format = "%Y-%m-%d %H:%M:%S")
pe_data$timestamp <- strptime(pe_data$timestamp, "%Y-%m-%d")
pe_data <- add_timeNumber(pe_data)

# add extra columns from the experimental design
pe_data <- add_exp_des_col(data = pe_data, d_exp_des = d_exp,
                           data_unit = "unit",
                           d_exp_unit = "new_unit",
                           col_add = c("rowNum", "colNum", "block", "cross"))

# PlotId is the combination of row and col information
pe_data$plotId <- paste0(paste0("c", pe_data$colNum), paste0("r", pe_data$rowNum))

# arrange the columns in a certain order
pe_data <- pe_data %>% select(timeNumber, timestamp, block, rowNum, colNum, plotId,
                              cross, genotype, Digital_biomass, Height,
                              Leaf_angle, Leaf_area, Leaf_area_index,
                              Leaf_area_projected, Leaf_inclination,
                              Light_penetration_depth) %>% rename(timePoint = timestamp)

# PE pipeline: Median calculation ----
pe_data <- median_computation(pe_data)

# PE pipeline: remove tp with too high missing values ----

plot_trend(data = pe_data, trait = "Height")
plot_trend(data = pe_data, trait = "Leaf_area")
plot_trend(data = pe_data, trait = "Digital_biomass")

prop_non_miss <- timepoint_prop_non_missing(pe_data)
tp_rem <- names(prop_non_miss[prop_non_miss < 0.3])
if(length(tp_rem) > 0){
  pe_data <- pe_data[!(pe_data$timePoint %in% tp_rem), ]
}

plot_trend(data = pe_data, trait = "Height")
plot_trend(data = pe_data, trait = "Leaf_area")
plot_trend(data = pe_data, trait = "Digital_biomass")

# PE pipeline: remove outliers ----

pe_data <- outlier_boxplot_detect(pe_data)

plot_trend(data = pe_data, trait = "Height")
plot_trend(data = pe_data, trait = "Leaf_area")
plot_trend(data = pe_data, trait = "Digital_biomass")

# PE pipeline: trim the time series ----

# for example remove days at the beginning or the end

sel_tp <- sort(unique(pe_data$timePoint)) # all tp in the TS
sel_tp <- sel_tp[-c(17:20)] # remove four last days

pe_data <- pe_data[pe_data$timePoint %in% sel_tp, ]

plot_trend(data = pe_data, trait = "Height")
plot_trend(data = pe_data, trait = "Leaf_area")
plot_trend(data = pe_data, trait = "Digital_biomass")

# PE pipeline: Creation of TP object ----

TP_PE <- createTimePoints(dat = pe_data,
                          experimentName = "Exp_51_Kenin_Keni_PE",
                          genotype = "genotype",
                          timePoint = "timePoint",
                          plotId = "plotId",
                          rowNum = "rowNum", colNum = "colNum")

# PE pipeline: Spatial adjustment ----

PE_adjusted <- spatial_adjustment(TP = TP_PE,
                                  traits = c("Height", "Leaf_area", "Digital_biomass"),
                                  timePoints = NULL,
                                  extraFixedFactors = "block",
                                  geno.decomp = NULL,
                                  what = "fixed",
                                  useCheck = FALSE,
                                  useRepId = FALSE,
                                  engine = "SpATS",
                                  spatial = FALSE,
                                  quiet = TRUE)

plot_PE_adj <- PE_adjusted$plot_res
geno_PE_adj <- PE_adjusted$geno_res
comp_mon <- PE_adjusted$comp_monitor
h2_PE <- PE_adjusted$h2_res

plot_trend(data = plot_PE_adj, trait = "Height_corr")
plot_trend(data = plot_PE_adj, trait = "Digital_biomass_corr")

plot_trend(data = geno_PE_adj, trait = "Leaf_area_pred", genotype = TRUE)
plot_trend(data = geno_PE_adj, trait = "Height_pred", genotype = TRUE)

# save the data
setwd(ref_wd)
save(plot_PE_adj, file = "PE_plot_adjusted.RData")
save(geno_PE_adj, file = "PE_geno_adjusted.RData")
save(h2_PE, file = "h2_PE.RData")

############################ Weather data comp #################################
# Weather data full TS computation ----
# get the weather data time series information
setwd(ref_wd)
load(file = "sensor_data_processed.RData")
load(file = "weather_data_processed.RData")

wth_data <- wth_data_proc(wth_data = wth_data, sensor_data = sensor_data)

# reduce the values over days
wth_data <- wth_data %>% group_by(dmy) %>% summarise(T_min = min(Temp, na.rm = TRUE),
                                                     T_max = max(Temp, na.rm = TRUE),
                                                     T_av = mean(Temp, na.rm = TRUE),
                                                     RH_av = mean(RH, na.rm = TRUE),
                                                     VPD_av = mean(VPD, na.rm = TRUE),
                                                     SR_av = mean(SR, na.rm = TRUE),
                                                     WS_av = mean(WS, na.rm = TRUE))
colnames(wth_data)[1] <- "timePoint"
wth_data$timePoint <- as.character(wth_data$timePoint)

# save wth data results
save(wth_data, file = "weather_data_full_TS.RData")

############################ merge LC + PE + wth ###############################
# load the adjusted data ----
setwd(ref_wd)
load(file = "LC_plot_adjusted.RData")
load(file = "LC_geno_adjusted.RData")

load(file = "PE_plot_adjusted.RData")
load(file = "PE_geno_adjusted.RData")

load(file = "weather_data_full_TS.RData")

# merge - plot data ----
plot_LC_adj <- plot_LC_adj[, -1]
plot_PE_adj <- plot_PE_adj[, -1]

# plot_data <- full_join(x = plot_LC_adj, y = plot_PE_adj, by = c("timePoint", "plotId"))
plot_data <- full_join(x = plot_LC_adj, y = plot_PE_adj,
                       by = intersect(colnames(plot_LC_adj), colnames(plot_PE_adj)))

# add wth data
plot_data$timePoint <- as.character(plot_data$timePoint)
plot_data <- left_join(x = plot_data, y = wth_data, by = "timePoint") %>%
  arrange(timePoint)

# save final data
save(plot_data, file = "plot_data.RData")

# merge - geno data ----
geno_LC_adj <- geno_LC_adj[, -1]
geno_PE_adj <- geno_PE_adj[, -1]

# geno_data <- full_join(x = geno_LC_adj, y = geno_PE_adj, by = c("timePoint", "genotype"))
geno_data <- full_join(x = geno_LC_adj, y = geno_PE_adj,
                       by = intersect(colnames(geno_LC_adj), colnames(geno_PE_adj)))

# add wth data
geno_data <- left_join(x = geno_data, y = wth_data, by = "timePoint") %>%
  arrange(timePoint)

# save final data
save(geno_data, file = "geno_data.RData")

# form a TP object for diagnostic plot ----
setwd(ref_wd)
load(file = "plot_data.RData")

TP <- createTimePoints(dat = plot_data,
                       experimentName = "Exp_51_Kenin_Keni",
                       genotype = "genotype",
                       timePoint = "timePoint",
                       plotId = "plotId",
                       rowNum = "rowId", colNum = "colId")

#### day 3 ----

tp <- 3

plot(TP, traits = "max_TR",
     plotType = "layout",
     timePoints = tp)

plot(TP, traits = "max_TR_corr",
     plotType = "layout",
     timePoints = tp)

plot(TP, traits = "total_auc",
     plotType = "layout",
     timePoints = tp)

plot(TP, traits = "total_auc_corr",
     plotType = "layout",
     timePoints = tp)

plot(TP, traits = "Height",
     plotType = "layout",
     timePoints = tp)

plot(TP, traits = "Height_corr",
     plotType = "layout",
     timePoints = tp)

plot(TP, traits = "Digital_biomass",
     plotType = "layout",
     timePoints = tp)

plot(TP, traits = "Digital_biomass_corr",
     plotType = "layout",
     timePoints = tp)

