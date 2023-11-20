library(readxl)

### Exp41-PM ###
load("./data/Exp22-Mulch_Cereals/m.lc.RData")
load("./data/Exp22-Mulch_Cereals/meta.d.RData")
load("./data/Exp22-Mulch_Cereals/t.rh.ws.RData")
load("./data/Exp22-Mulch_Cereals/solRAD.RData")
load("./data/Exp22-Mulch_Cereals/pe.df.RData")
sensor.unit.df <- read.csv("./data/Exp22-Mulch_Cereals/sensor_unit_map.csv")

HTP_data_Exp22_Cereals <- list(m.lc = m.lc, meta.d = meta.d, 
                               clm.df = as.data.frame(rbind(t.rh.ws, solRAD)),
                               sensor.unit.df = sensor.unit.df, 
                               pe.df = pe.df)
save(HTP_data_Exp22_Cereals, 
     file = "./data/Exp22-Mulch_Cereals/HTP_data_Exp22_Cereals.RData")
