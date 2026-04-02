#Boxplot

DefaultAssay(onlyBE2) <- "RNA"
Idents(object = onlyBE2) <- "ARQExp"
boxdata = FetchData(onlyBE2, c("ARQExp", "ARQ", "Igf1r", "Fos", "Jak2", "Mapk13"))
tail(boxdata,6)

tiff(file = "onlyBE2 hARtg Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Mapk13, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()

#Volinplot with Mean +- s.d.(standard deviation)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

DefaultAssay(TriplevDouble.combined1.BELEMETPos) <- "RNA"
Idents(object = TriplevDouble.combined1.BELEMETPos) <- "stim"

tiff(file = "TriplevDouble.combined1.BELEMETPos Xpo1 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Xpo1", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 1) #you can change/delete y.max = 
+ NoLegend() 
+stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
