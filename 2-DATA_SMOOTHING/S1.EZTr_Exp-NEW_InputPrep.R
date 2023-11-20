setwd(dir = "C:/Users/2021lg003/Documents/Sorghum/Experiments/EXP_INDE_2023/ANALYSIS_ALL_DATASET/DATA_SMOOTHING")
load("./HTP_data.RData")
m.lc <- HTP_data$m.lc
meta.d <- HTP_data$meta.d
t.rh.ws <- HTP_data$t.rh.ws
solRAD <- HTP_data$solRAD
pe.df <- HTP_data$pe.df

library(dplyr)
library(lubridate)

### Create meta.d ###
allsec.d <- read.csv("./data/PERIODE2-01-10-10-10/meta.d_raw.csv", sep = ";")
meta.d<-allsec.d

write.csv(meta.d, "data/PERIODE2-01-10-10-10/meta.d.csv")
save(meta.d, file = "data/PERIODE2-01-10-10-10/meta.d.RData")

### Create m.lc ###
# cereals original LC data #
all.lc <- read.csv("data/PERIODE2-01-10-10-10/m.lc_raw_PERIODE2.csv", sep = ";")
new.lc<-as.data.frame(matrix(nrow = nrow(all.lc), ncol = ncol(m.lc)))
colnames(new.lc) <- colnames(m.lc)
#new.lc$unit <- paste0(all.lc$Row,"-1",":",all.lc$Column)
new.lc$unit <- all.lc$Row
new.lc$timestamp <- all.lc$Timestamp     ###ici pb pour appeler all.lc$Timestamp  si on ne met pas de separateur en amont, 
new.lc$timestamp <- dmy_hm(new.lc$timestamp) #attention au format de cellule sous Excel AAAA.MM.JJ HH:MM
new.lc$Mass..g.<-all.lc$Weight_g
for (i in 1:nrow(new.lc)) { 
  r.ind <- which(meta.d$unit==new.lc$unit[i])
  new.lc$genotype[i]<-unique(meta.d$Genotype[r.ind])
  new.lc$g_alias[i]<-unique(meta.d$G..Alias[r.ind])
  new.lc$treatment[i]<-unique(meta.d$Treatment[r.ind])
}
m.lc <- new.lc
write.csv(m.lc, "data/PERIODE2-01-10-10-10/m.lc_PERIODE2.csv")
save(m.lc, file = "data/PERIODE2-01-10-10-10/m.lc_PERIODE2.RData")



### Create Climate data ###

clm <- read.csv("data/PERIODE2-01-10-10-10/climate_raw_PERIODE2.csv", sep = ";") ##ici, pb de colonne : clm.n <- clm[,1:7] si on ne met pas de separateur ;
clm.n <- clm[,1:7]

clm.n$Date <- dmy(clm.n$Date)
clm.n$RH <- substr(clm.n$RH, 1, 2)  
clm.n$RH <- as.numeric(clm.n$RH)  

write.csv(clm.n, "data/PERIODE2-01-10-10-10/climate_PERIODE2.csv")
save(clm.n, file = "data/PERIODE2-01-10-10-10/climate_PERIODE2.RData")
  
### Exp-Layla2a ###
IRD_data_Exp_NEW <- list(m.lc = m.lc, meta.d = meta.d, climate = clm.n)

save(IRD_data_Exp_NEW, file = "data/PERIODE2-01-10-10-10/IRD_data_Exp_NEW.RData")

###########################################################



### Create t.rh.ws and solRAD ###
# read original weather files of cereals mulch exp#
# t.rh.d <- read.csv("data/Exp22-Mulch_Cereals/Etr-LAI-cereals1_org_data_climate_T-RH.csv")
# solrad.ws.d <- read.csv("data/Exp22-Mulch_Cereals/Etr-LAI-cereals1_org_data_climate_Batt_min-Rain-SlrW-Wind.csv")
# 
# new.t <- t.rh.d[,c(1,2,4)]
# new.rh <- t.rh.d[,c(1,2,3)]
# new.ws <- solrad.ws.d[,c(1,2,6)]
# new.solrad <- solrad.ws.d[,c(1,2,5)]
# 
# colnames(new.ws)<-names(t.rh.ws)[c(1,3,4)]
# colnames(new.solrad) <- colnames(new.ws)
# colnames(new.t)<-colnames(new.rh)<-colnames(new.ws)
# 
# new.t$variable <- rep("Temperature (Â°C)", nrow(new.t))
# new.rh$variable <- rep("Relative humidity (%)", nrow(new.rh))
# new.ws$variable <- rep("Windspeed average (m/s)", nrow(new.ws))
# new.solrad$variable <- rep("Solar radiation (W/(s*mÂ²))", nrow(new.solrad))
# 
# new.t <- new.t[,c(1,4,2,3)]
# new.rh <- new.rh[,c(1,4,2,3)]
# new.ws <- new.ws[,c(1,4,2,3)]
# new.solrad <- new.solrad[,c(1,4,2,3)]
# 
# new.t.rh.ws <- as.data.frame(rbind(new.t, new.rh))
# newsolRAD <- as.data.frame(rbind(new.ws, new.solrad))
# 
# t.rh.ws <- new.t.rh.ws
# solRAD <- newsolRAD
# 
# write.csv(t.rh.ws,
#           "D:/MULCH_TRTMNT DATA/CEREALS_RESULTS/Mulch_Cereals/data/InputFiles/t.rh.ws.csv")
# save(t.rh.ws,
#      file = "D:/MULCH_TRTMNT DATA/CEREALS_RESULTS/Mulch_Cereals/data/InputFiles/t.rh.ws.RData")
# 
# write.csv(solRAD,
#           "D:/MULCH_TRTMNT DATA/CEREALS_RESULTS/Mulch_Cereals/data/InputFiles/solRAD.csv")
# save(solRAD,
#      file = "D:/MULCH_TRTMNT DATA/CEREALS_RESULTS/Mulch_Cereals/data/InputFiles/solRAD.RData")
# 
# ### create pe.df ###
# new.allsecPE <- allsec.d
# substr(new.allsecPE$BCFG.data_Sector, 4, 4)<-":"
# names(new.allsecPE)[c(1,6)]<-c("Sector","G..Alias")
# new.allsecPE <- new.allsecPE[,c(1:8,13,9,14,10,11,15,12,16)]
# new.allsecPE$Treatment<-substr(new.allsecPE$Treatment,5,nchar(new.allsecPE$Treatment))
# new.allsecPE$Genotype <- factor(new.allsecPE$Genotype)
# new.allsecPE$G..Alias <- factor(new.allsecPE$G..Alias)
# new.allsecPE$Replicates <- factor(new.allsecPE$Replicates)
# 
# pe.df <- new.allsecPE
# 
# write.csv(pe.df,
#           "D:/MULCH_TRTMNT DATA/CEREALS_RESULTS/Mulch_Cereals/data/InputFiles/pe.df.csv")
# save(pe.df,
#      file = "D:/MULCH_TRTMNT DATA/CEREALS_RESULTS/Mulch_Cereals/data/InputFiles/pe.df.RData")
# 
# HTP_data_Exp22_Cereals <- list(m.lc=m.lc, pe.df=pe.df, meta.d=meta.d,
#                                solRAD = solRAD, t.rh.ws=t.rh.ws)
# save(HTP_data_Exp22_Cereals, 
#      file = "D:/MULCH_TRTMNT DATA/CEREALS_RESULTS/Mulch_Cereals/data/InputFiles/HTP_data_Exp22_Cereals.RData")
# load("D:/MULCH_TRTMNT DATA/CEREALS_RESULTS/Mulch_Cereals/data/InputFiles/HTP_data_Exp22_Cereals.RData")
# allData <- HTP_data_Exp22_Cereals
