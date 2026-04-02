m <- read.table("../codes_omics_data/peak calling/intermediate_data/1_merge.txt")
candidate25 <- read.table("../codes_omics_data/peak calling/intermediate_data/25peaks.txt")
pdf(file =  '25peaks.pdf', height = 12, width = 12)
par(mfrow=c(5,5))
par(mar=c(4.5, 4.5, 2, 1))
for (i in 1:25) {
  s <- candidate25[i,1]
  e <- candidate25[i,2]
  s1 <- s-2000
  e1 <- e+2000
  subm <- m[m[,1]<=e1 & m[,2]>=s1,]
  subm[1,1] <-s1
  subm[nrow(subm),2] <- e1
  x <- s1:e1
  y1 <- rep(0, e1-s1)
  y2 <- rep(0, e1-s1)
  
  for (j in 1:nrow(subm)) {
    y1[(subm[j,1]-s1+1):(subm[j,2]-s1+1)]=subm[j,3]
    y2[(subm[j,1]-s1+1):(subm[j,2]-s1+1)]=subm[j,4]
  }
  plot(x,y1,'l',xlab = i, ylab = "Read coverage", cex.lab = 1.6, cex.axis = 1.4)
  lines(x,y2,'l', col = 'blue')
  abline(v=s, lty= 2)
  abline(v=e, lty= 2)
}
dev.off()

