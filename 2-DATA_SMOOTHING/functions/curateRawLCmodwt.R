
curateRawLC.Et <- function(x, y) {
  
  x = LC.MAT.f; y = meta.LCDF
  
  base.d <- x[,-1]
  
  data.RmOut <- base.d 


# stat.mat <- as.data.frame(apply(base.d, 2, function(x) {quantile(x, na.rm = TRUE)}))
  allStats <- as.data.frame(matrix(NA, nrow = 11, ncol = ncol(base.d)))
  rownames(allStats) <- c("mean", "sd", "median", "trimmed", "mad", "min", 
                          "max", "range", "skew", "kurtosis", "se")
  colnames(allStats) <- colnames(base.d)
  
  for(i in 1:ncol(base.d))
  {
    tmp <- c(describe(base.d[ ,i])[-c(1:2)])
    allStats[,i] <- c(unlist(describe(base.d[ ,i])[-c(1:2)]))
  }
  
  
  # Function to get lower and upper limits
  maxL <- function(x) {(mean(as.numeric(x))+1.5*sd(as.numeric(x)))}
  minL <- function(x) {(mean(as.numeric(x))-1.5*sd(as.numeric(x)))}
  
  # Sector-IDs with values beyond limits
  min.OL0<- which(allStats[6, ] < minL(allStats[6, ]))
  max.OL100 <- which(allStats[7, ] > maxL(allStats[7, ]))
  max.OLskw <- which(allStats[9, ] > maxL(allStats[9, ]))
  min.OLskw <- which(allStats[9, ] < minL(allStats[9, ]))
  max.OLkurt <- which(allStats[10, ] > maxL(allStats[10, ]))
  min.OLkurt <- which(allStats[10, ] < minL(allStats[10, ]))
  max.OLse <- which(allStats[11, ] > maxL(allStats[11, ]))
  min.OLse <- which(allStats[11, ] < minL(allStats[11, ]))
  
  olSECs <- unique(c(min.OL0, max.OL100, max.OLskw, min.OLskw, 
                     max.OLkurt, min.OLkurt, max.OLse, min.OLse))
  # olSECnms <- as.character(LC.MAT.raw$old_unit[olSECs])

# Remove outliers from data  
base.d <- base.d[ ,-olSECs]
data.RmOut <- data.RmOut[ ,-olSECs]
y.tmp <- y[y$unit %in% colnames(base.d), ]

for ( i in 1:ncol(base.d)) {
  data <- as.numeric(c(base.d[, i]))
  data[data < 0] <- NA # if less than 0, it is wrong.
  res.dwt <- modwt(data, filter = "haar", n.level = 3, boundary = "periodic", fast = FALSE)
  bp.s1 <- boxplot.stats(res.dwt@W$W1); tf1 <- res.dwt@W$W1 %in% bp.s1$out#; tf1 <- MyFun(tf1)
  bp.s2 <- boxplot.stats(res.dwt@W$W2); tf2 <- res.dwt@W$W2 %in% bp.s2$out#; tf2 <- MyFun(tf2)
  bp.s3 <- boxplot.stats(res.dwt@W$W3); tf3 <- res.dwt@W$W3 %in% bp.s3$out#; tf3 <- MyFun(tf3)
  tf <- tf1 | tf2 | tf3
  data.RmOut[tf, i] <- NA
  
  # jpeg(paste0(opPATH, "Plot_dwt.OL/", y.tmp$old_unit[i], ".jpeg"))
  # plot(data, cex = 1, lwd = 0.8)
  # points(x = (1:length(data))[!tf], y = data[!tf], col = "red", pch = 20, cex = 0.8)
  # dev.off()  
}
# write.csv(data.RmOut, file = "./RE-ANALYSIS_29.06/S1_2 results/allSECs.dwt.OL_s1.csv")

rmOut <- as.data.frame(t(data.RmOut))

rmOut.DF <- as.data.frame(cbind(y.tmp, rmOut))
colnames(rmOut.DF) <- c(colnames(y.tmp), as.character(x$TS))

####### Imputation #######
data.Imp <- data.frame(apply(data.RmOut, 2, na.approx, na.rm = F))

data.Imp.df <- as.data.frame(t(data.Imp))

imputed.DF <- as.data.frame(cbind(y.tmp, data.Imp.df))
colnames(imputed.DF) <- c(colnames(y.tmp), as.character(x$TS))

# write.csv(imputed.DF, paste0(opPATH, "LC.MAT.BOX.OL_s1_IMPUTED.csv"))

# # There still could remain columns with maximum missing, hence run the below script
# na.list <- as.numeric(apply(imputed.DF[,-c(1:5)], 1, 
#                             FUN = function(x) {sum(is.na(x))}))
na.list <- c()
for (c in 1:nrow(imputed.DF)) {
  cnt <- sum(is.na(imputed.DF[c, 6:ncol(imputed.DF)]))
  na.list<-c(na.list, cnt)
}

# keep sectors which have more than 30% of values
na.list.G.Locs <- which(na.list > ceiling(0.3*dim(imputed.DF[,-c(1:5)])[2])) 

# replace 'na.list.G.Locs' only with 0 if NA so that ETr can be extracted
if(length(na.list.G.Locs) > 0)
{
  imputed.DF.tmp <- imputed.DF
  imputed.DF.tmp <- imputed.DF.tmp[-(na.list.G.Locs+5), ]
} else {imputed.DF.tmp <- imputed.DF}

imputed.DF.final <- imputed.DF.tmp

interp.ip <- imputed.DF.final[ ,6:ncol(imputed.DF.final)]

interp.df <- as.data.frame(t(interp.ip))

interp.df.op <- interp.df

for (c in 1:nrow(interp.df)) {
  interp.df.op <- na.aggregate.default(interp.df[c, ])
  }

# interp.df.op <- as.data.frame(apply(interp.df, 1, na.aggregate.default))

imputed.DF.final[ ,6:ncol(imputed.DF.final)] <- interp.df.op

return(imputed.DF.final)
}













# # Remove outlier again
# data.RmOut.2 <- as.data.frame(data.RmOut)
# for (i in 1:ncol(base.d) ) {
#   data <- data.RmOut[, i]
#   # plot(data)
#   # points(data.im, col = "red", pch = 20, cex = .8)
#   data.im <- na.locf(data, na.rm = FALSE)
#   DIF <- diff(data.im)
#   Threshold <- quantile(DIF, prob = 1 - (4 / length(DIF)), na.rm = T)
#   #sum(DIF > Threshold, na.rm = T)
#   tf <- DIF > Threshold
#   tf[is.na(tf)] <- FALSE # NA should be FALSE
#   tf <- c(FALSE, tf) # NA should be FALSE
#   num <- (1:length(tf))[tf]
#   num.rm <- as.numeric(c(sapply(num, "+", 0:12)))
#   data.RmOut.2[num.rm, i] <- NA
#   # 
#   # jpeg(paste0(opPATH, "Plot_EMPRCL.OL/", y.tmp$old_unit[i], ".jpeg"))
#   # plot(data)
#   # points(num.rm, data[num.rm], col = "red", pch = 20, cex = 0.9)
#   # dev.off() 
#   
#   Sys.sleep(0.3)
# }
# 
# ext <- which(!rownames(data.RmOut.2) %in% rownames(data.RmOut))
# data.RmOut.2 <- data.RmOut.2[-ext,]
# 
# 
# # Imputation
# data.Imp <- data.frame(apply(data.RmOut.2, 2, na.approx, na.rm = F))
# rownames(data.Imp) <- rownames(data.RmOut.2)
# data.final <- as.data.frame(apply(data.Imp, 2, na.aggregate))
# rownames(data.final) <- rownames(data.RmOut.2)
# 
# # ET values
# data.et <- data.final
# data.et[1, ]<-0
# 
# for(i in 1:ncol(data.final))
# {
#   data.et[-1,i] <- -diff(data.final[, i])
#   
#   # jpeg(paste0(opPATH, "Plot_ET.final/", y.tmp$old_unit[i], ".jpeg"))
#   # plot(-diff(data.final[, i]), type = "p")
#   # points(-diff(data.final[, i]), type = "l")
#   # dev.off() 
# }
# 
# rownames(data.et) <- x$TS
# etMeta <- as.data.frame(cbind(y.tmp, t(data.et)))
# 
# return(etMeta)

# }
