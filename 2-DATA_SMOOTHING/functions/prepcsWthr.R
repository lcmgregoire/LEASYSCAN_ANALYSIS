
prepcsWthr <- function(x, y) {
  
  
  temp.mat <- matrix(x[ ,y], nrow = (24*60), 
                     ncol = nrow(x)/(24*60), byrow = FALSE)
  
  # remove outliers from each date
  for(i in 1:ncol(temp.mat))
  {
    df <- temp.mat[,i]
    
    # if (y==2|y==3) # if wthr var is Temp/Rel.Hum, replace 0 with NA
    # {df[df == 0.0] <- NA}
    
    temp.mat[,i] <- df
    
    ol <- boxplot(temp.mat[,i], plot = FALSE)$out
    
    temp.mat[temp.mat[,i] %in% ol, i] <- NA
  }
  
  ########## MICE imputation ##########
  # mice_plot <- aggr(temp.mat, col=c('navyblue','yellow'),
  #                   numbers=TRUE, sortVars=FALSE,
  #                   labels=names(temp.mat), cex.axis=.7,
  #                   gap=3, ylab=NULL)
  
  # Impute the missing values.
  # imputed_Data <- mice(temp.mat, m=5, maxit = 10, method = 'pmm', 
  #                      seed = 500, printFlag = FALSE)
  # 
  # # check imputed values
  # OP <- mice ::complete(imputed_Data, 5, include = FALSE)
  # 
  # OP.smth <- OP
  # 
  # for (s in 1:ncol(OP)) {
  #   opdf <- as.data.frame(OP[ ,s])
  #   
  #   if(sum(is.na(opdf)) >= 0.9*nrow(opdf)){
  #     OP.smth[ ,s] <- opdf[,1]
  #   }else{
  #     smth.opdf <- smooth.spline(opdf)
  #     OP.smth[ ,s] <- smth.opdf$y
  #   }
  #   
  #   }
  # 
  # OP.vec <- as.vector(as.matrix(OP.smth))
  # OP.vec <- na.aggregate.default(OP.vec)
  
  ############################################################
  
  OP.vec <- as.vector(as.matrix(temp.mat))

  OP.vec <- na.aggregate.default(OP.vec)
  
  return(OP.vec)
}