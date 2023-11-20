filterLCExtremeCols <- function(x, y) {
  
  des.stats.ip <- as.data.frame(t(x[,-c(1:6)]))
  colnames(des.stats.ip) <- x$unit
  
  ol.cnts <- matrix(nrow = 1, ncol = ncol(des.stats.ip))
  colnames(ol.cnts) <- x$unit
  
  for (i in 1:ncol(des.stats.ip)) {
    des.stats.ip[,i][is.nan(des.stats.ip[,i])]<-NA 
    ol.cnts[i]<-length(boxplot(des.stats.ip[,i], plot = FALSE)$out)
  }
  
  # des.stats <- sapply(na.omit(des.stats.ip),
  #                     FUN = function(x){length(boxplot(des.stats.ip, plot = FALSE)$out)})
  # 
  # q0 <- which(des.stats[1, ] %in% boxplot(des.stats[1, ], plot = FALSE)$out)
  # q25 <- which(des.stats[2, ] %in% boxplot(des.stats[2, ], plot = FALSE)$out)
  # q50 <- which(des.stats[3, ] %in% boxplot(des.stats[3, ], plot = FALSE)$out)
  # q75 <- which(des.stats[4, ] %in% boxplot(des.stats[4, ], plot = FALSE)$out)
  # q100 <- which(des.stats[5, ] %in% boxplot(des.stats[5, ], plot = FALSE)$out)
  
  # err.cols <- intersect(intersect(intersect(intersect(q25, q50), q75), q100), q0)
  
  err.cols <- which(ol.cnts[1, ] > 0.5*dim(des.stats.ip)[1])
  
  err.sec.nm <- colnames(ol.cnts)[err.cols]
  
  err.sec.meta <- y[y$unit %in% err.sec.nm, ]
  
  list(err.sec.META = err.sec.meta, err.sec.NM = err.sec.nm)
}