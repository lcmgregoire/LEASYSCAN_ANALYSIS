
wthr.vars <- read.csv("results/Sorghum_weather.DF.agg15min_MACRO vals.csv")
wthr.vars <- wthr.vars[,-1]                      

SVP <- 610.7*(10^(7.5*wthr.vars[ ,2]/(237.3+wthr.vars[ ,2])))
VPD <- ((1 - (wthr.vars[ ,3]/100))*SVP)/1000
wthr.vars[ ,4] <- VPD

write.csv(wthr.vars, "results/Sorghum_Wthr_agg15min_MACROvals_VPD.csv")
