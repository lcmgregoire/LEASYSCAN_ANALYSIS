calculateETref <- function(x) {
  
  wthr.df <- as.data.frame(x)
  
  ## Check date and number format in wthr.df and then use the functions
  wthr.df$TS<-ymd_hms(wthr.df$TS) 
  
  wthr.df[ ,2:ncol(wthr.df)] <- apply(wthr.df[ ,2:ncol(wthr.df)], 2, 
                                      function(x) {as.numeric(as.character(x))})
  
  ## extract observations for EACH DATE 
  
  wthr.df$date <- date(wthr.df$TS)
  
  ## extract observations for EACH TIME 
  
  wthr.df$time <- strftime(wthr.df$TS, format="%H:%M:%S", tz="UTC")
  
  # wthr.df %>% group_by(time) %>% mutate(Tmax = summarise(Value = max(wthr.df$T))) 
  
  ## Group by date
  wthr.df1<-wthr.df %>%
    arrange(date, time) %>%
    group_by(date) %>%
    mutate(Tmax = max(Temp), Tmin = min(Temp)) 
  
  
  ## Calculate the Penman-Monteith Equation
  
  ########################## ESPACE DE TRAVAIL
  ####### le calcul d'ETREF est faux. REVOYONS L'EQUATION DE PENMAN
  ##voir le fichier ETREF MANIP SEAUX pour voir les formules correctement + FAO : https://www.fao.org/3/X0490E/x0490e06.htm#chapter%202%20%20%20fao%20penman%20monteith%20equation
  
  
  #ET = (0.408 * ∆ * (Rn - G) + γ (900/(T+273)) * u2 * (es - ea)) /  ∆  + γ (1+0.34u2)
  
 #where
  
# ETo reference evapotranspiration [mm day-1],
 #Rn net radiation at the crop surface [MJ m-2 day-1],  ==== The net radiation, Rn, is the difference between incoming and outgoing radiation of both short and long wavelengths. It is the balance between the energy absorbed, reflected and emitted by the earth's surface or the difference between the incoming net shortwave (Rns) and the net outgoing longwave (Rnl) radiation (Figure 15). Rn is normally positive during the daytime and negative during the nighttime. The total daily value for Rn is almost always positive over a period of 24 hours, except in extreme conditions at high latitudes. 
 #G soil heat flux density [MJ m-2 day-1],              ==== In making estimates of evapotranspiration, all terms of the energy balance (Equation 1) should be considered. The soil heat flux, G, is the energy that is utilized in heating the soil. G is positive when the soil is warming and negative when the soil is cooling. Although the soil heat flux is small compared to Rn and may often be ignored, the amount of energy gained or lost by the soil in this process should theoretically be subtracted or added to Rn when estimating evapotranspiration.
# T mean daily air temperature at 2 m height [°C],
 #u2 wind speed at 2 m height [m s-1],                  ==== ws
 #es saturation vapour pressure [kPa],                  ==== es = 0.6108*exp(17.27*t/(t + 273.3))
 #ea actual vapour pressure [kPa],                      ==== ea <- rh/100*es
# es - ea saturation vapour pressure deficit [kPa],
 #∆ slope vapour pressure curve [kPa °C-1],             ==== 4098(0.6108*exp(17.27*t/(t+237.3))))/(t+273.3)^2 
# gamma psychrometric constant [kPa °C-1].              ====
  
  ####################################################
  ### RETOUR AU PIPELINE
  
  # del <- (4098(0.6108*exp(17.27*t/(t+237.3))))/(t+273.3)^2 ###
  # 
  # gammaa <- 0.63189 # P at 500m msl (meters sea level?) = 95.52 kPa & gamma = i.e.0.000665*95.52  
  # 
  # es = 0.6108*exp(17.27*t/(t + 273.3))
  # 
  # ea <- rh/100*es
  # 
  # Rng <- ((slr.rad*15*60/1000000)*2.5)
  #  
  # ETref <- ((0.408*del*Rng) + 0.063189*((900/96)*ws*(es-ea))/(t+273))/(del+0.063189*(1+0.34*ws))
  
  wthr.df1$ETref <- NA
  
  for(i in 1:nrow(wthr.df1))  {
    t <- wthr.df1$Temp[i]
    ws <- wthr.df1$WS[i]+0.000001 ### pourquoi ajouter 0.00001? (unité pour Penman : m/s)
    tmax <- wthr.df1$Tmax[i]
    tmin <- wthr.df1$Tmin[i]
    rh <- wthr.df1$RH[i]
    slr.rad <- wthr.df1$SR[i]*10^-6 ### pourquoi multiplier par *10^6 ? (unité pour Penman : MJ/s * m²) Rn-G ## kar 2020 dit qu'on a utiliser les stations Campbell, donc en W/m²?
    # si on a des W =  J/s --- on veut ajouter MJ --1MJ = 10e-6 MJ
    del <- (4098*(0.6108*exp(17.27*t/(t+237.3))))/(t+273.3)^2
    
    gammaa <- 0.63189  # P at 500m msl = 95.52 kPa
    
    es = 0.6108*exp(17.27*t/(t + 273.3))
    
    ea <- rh/100*es
    
    Rng <- ((slr.rad*15*60)*2.5) ##15*60 pour convertir en /s : *15 min * 60 sec, 
    ## pourquoi /1 000 000?  pourquoi *2.5? 
    
    ETref <- ((0.408*del*Rng) + 0.063189*((900/96)*ws*(es-ea))/(t+273))/(del+0.063189*(1+0.34*ws)) ## 96 is the constant/
   #ETref <- ((0.408*∆*Rng) + gamma*((900/96)*ws*(es-ea))/(t+273))/(del+0.063189*(1+0.34*ws))
    if(ETref == 0.0){
      
      wthr.df1$ETref[i] <- NA
      
    } else {wthr.df1$ETref[i] <- ETref}
  }
  
  wthr.df1$ETref <- na.approx(wthr.df1$ETref)
  
  wthr.df1$ETref<-round(wthr.df1$ETref, 3)
  
  return(wthr.df1)
}
