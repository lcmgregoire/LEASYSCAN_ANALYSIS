smoothETr <- function(x) {
  
  etrOP <- etrIP <- x
  
  n <- ncol(etrIP)
  
 for(i in 1:nrow(etrIP))
  {
   dat = as.data.frame(matrix(nrow = n, ncol = 2))
   names(dat) <- c("y", "x")
   dat$y = t(etrIP[i,]) 
   dat$x = 1:n
      
   outliers <- which(dat$y < -0.05)
   dat$y[outliers] = NA
   if(is.na(dat$y[1])){dat$y[1]==0.0}
   if(is.na(dat$y[n])){dat$y[1]==0.0}
   dat$y <- na.spline(dat$y)
   
   loessData <- data.frame(
      x = 1:n,
      yhat = predict(loess(y~x, dat, span = 0.10)),
      method = "loess()" )
   
   # loessData$y <- dat$y
   
   etrOP[i, ]<-loessData$yhat
  }
  return(etrOP)
  
  # e21op<-x
  # 
  # e21ip<-(e21op)*1000 ## was not working for values ranged between 0-1, hence multiplied with 10^3
  # 
  # e21ip[e21ip < -100]<--1 ## replace extreme -ve values (irrigation drainage) with -1
  # 
  # for(i in 1:nrow(e21ip))
  # {
  #   nlTS<- as.data.frame(nonLinearNoiseReduction(e21ip[i, ], 3, 12))
  #   names(nlTS)<-"x"
  #   
  #   nlTS$x <- c(e21ip[i, ])
  #   
  #   # Use O/P of 1 as I/P to Spline Smoothing
  #   
  #   spldf<-smooth.spline(nlTS)
  #   
  #   e21op[i, ]<-spldf$y
  #   
  # }
  # 
  # # plot(x=1:ncol(e21ip), y=(e21ip[10, ])/1000, col="black", type="p", main = "An example plot")
  # # lines(x=1:ncol(e21ip), y=(e21op[10, ])/1000, col="red")
  #  
  # return(e21op)
  
}