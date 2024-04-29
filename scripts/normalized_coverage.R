library(ggplot2)
library(data.table)


PerRefPATH = snakemake@input[[1]]
PerWindows = snakemake@input[[2]]

sampl = snakemake@params[[1]]
sp0 = snakemake@params[[2]]
outBase = snakemake@params[[3]]

SIZE = snakemake@config[["size"]]
STEP = snakemake@config[["step"]]
quality = snakemake@config[["quality"]]



depth_avg <- fread(PerWindows)
genome_Stats <- fread(PerRefPATH)
depth_avg$Ref <- sapply(strsplit(as.character(depth_avg$Chr), "_"), `[`, 1)
depth_avg$chr <- sapply(strsplit(as.character(depth_avg$Chr), "_"), `[`, 5)
depth_avg$Chr <- paste(depth_avg$Ref,depth_avg$chr,sep="_chr")

sp = strsplit(sp0,"_")[[1]][1]

rain.name = paste(outBase,sampl,"_",sp0,"_",quality,"_",SIZE,"by",STEP,"_",sp,"_normalized.png",sep="")
rain.name.m = paste(outBase,sampl,"_",sp0,"_",quality,"_",SIZE,"by",STEP,"_",sp,"_mt_normalized.png",sep="")
rain.name.2 = paste(outBase,sampl,"_",sp0,"_",quality,"_",SIZE,"by",STEP,"_",sp,"_2m_normalized.png",sep="")
f.plt = paste(outBase,sampl,"_",sp0,"_",quality,"_",SIZE,"by",STEP,"_normalized_windows.csv",sep="")
f.cont = paste(outBase,sampl,"_",sp0,"_",quality,"_",SIZE,"by",STEP,"_normalized_contigs.csv",sep="")

mdn <-genome_Stats[1,"Mean"]

perc5 = depth_avg[chr!="mt" & chr!="2m"]
perc5[,coverage.mean := format(round(coverage.mean/as.double(mdn), 4),nsmall=4)]

p <- ggplot(perc5, aes(x=window.end, y=as.double(coverage.mean),fill = chr,colour= chr)) +
    geom_point(shape = 20, size = 0.7,alpha=0.5) + scale_x_continuous(name="Genomic Position (bp)",
    labels = scales::scientific) + scale_y_continuous(name="Normalized Coverage Depth",
    limits=c(0, 5), breaks = seq(0,5,0.5)) + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_grid(. ~ chr, space="free_x", scales = "free_x",switch="x")+
    theme(panel.spacing.x = grid::unit(0, "cm")) + theme(text = element_text(size=10),
    axis.text.x = element_blank(),axis.ticks.x=element_blank()) +
    theme(panel.border = element_rect(colour = "grey",size = 0.1),strip.background = element_blank(),
    strip.text = element_text(angle =90),strip.placement = "outside") +
    geom_hline(aes(yintercept = 0.5)) + geom_hline(aes(yintercept = 1)) +
    geom_hline(aes(yintercept = 1.5)) + geom_hline(aes(yintercept = 2)) +
    geom_hline(aes(yintercept = 2.5)) + geom_hline(aes(yintercept = 3)) +
    geom_hline(aes(yintercept = 3.5)) + geom_hline(aes(yintercept = 4)) +
    ggtitle(paste(sampl,"_",sp," (",mdn,")",sep="")) +
    guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(rain.name,plot = p, width = 18, height = 12, units = "cm",bg = "white", dpi = 300)


perc5.m <- depth_avg[chr=="mt"]
perc5.m[,coverage.mean := format(round(coverage.mean/as.double(mdn), 4),nsmall=4)]
upr <- as.numeric(ceiling((as.double(max(perc5.m$coverage.mean)) + (genome_Stats[as.logical(genome_Stats[,c("Ref")]==sp),c("Max")]/mdn))/2))
p.m <- ggplot(perc5.m, aes(x=window.end, y=as.double(coverage.mean),fill = chr,colour= chr)) +
    geom_point(shape = 20, size = 0.7,alpha=0.5) + scale_x_continuous(name="Genomic Position (bp)",
    labels = scales::scientific) + scale_y_continuous(name="Normalized Coverage Depth",
    limits=c(0, upr), breaks = seq(0,upr,2)) + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_grid(. ~ chr, space="free_x", scales = "free_x",switch="x")+
    theme(panel.spacing.x = grid::unit(0, "cm")) + theme(text = element_text(size=10),
    axis.text.x = element_blank(),axis.ticks.x=element_blank()) +
    theme(panel.border = element_rect(colour = "grey",size = 0.1),strip.background = element_blank(),
    strip.text = element_text(angle =90),strip.placement = "outside") +
    ggtitle(paste(sampl,"_",sp," (",mdn,")",sep="")) +
    guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(rain.name.m,plot = p.m, width = 18, height = 12, units = "cm",bg = "white", dpi = 300)

perc5.2 <- depth_avg[chr=="2m"]
perc5.2[,coverage.mean := format(round(coverage.mean/as.double(mdn), 4),nsmall=4)]
upr2 <- as.double(ceiling((as.double(max(perc5.2$coverage.mean)) + (genome_Stats[as.logical(genome_Stats[,c("Ref")]==sp),c("Max")]/mdn))/2))
p.2 <- ggplot(perc5.2, aes(x=window.end, y=as.double(coverage.mean),fill = chr,colour= chr)) +
    geom_point(shape = 20, size = 0.7,alpha=0.5) + scale_x_continuous(name="Genomic Position (bp)",
    labels = scales::scientific) + scale_y_continuous(name="Normalized Coverage Depth",
    limits=c(0, upr), breaks = seq(0,upr,2)) + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_grid(. ~ chr, space="free_x", scales = "free_x",switch="x")+
    theme(panel.spacing.x = grid::unit(0, "cm")) + theme(text = element_text(size=10),
    axis.text.x = element_blank(),axis.ticks.x=element_blank()) +
    theme(panel.border = element_rect(colour = "grey",size = 0.1),strip.background = element_blank(),
    strip.text = element_text(angle =90),strip.placement = "outside") +
    ggtitle(paste(sampl,"_",sp," (",mdn,")",sep="")) +
    guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(rain.name.2,plot = p.2, width = 18, height = 12, units = "cm",bg = "white", dpi = 300)


perc5 <- merge(perc5,perc5.m,all.y=T,all.x=T)
perc5 <- merge(perc5,perc5.2,all.y=T,all.x=T)
write.csv(perc5[,-c(5,6)],f.plt,row.names=F)
perc.cont <- perc5[,.(Mean = format(round(mean(as.double(coverage.mean)),4),nsmall=4), SD = format(round(sd(as.double(coverage.mean)),4),nsmall=4),
  Min = min(as.double(coverage.mean)), Max = max(as.double(coverage.mean)),
  l_quantile = quantile(as.double(coverage.mean),probs=0.25),Median = median(as.double(coverage.mean)),
  u_quantile = quantile(as.double(coverage.mean),probs=0.75)),.(Chr)]
for(i in 1:length(perc.cont$Chr)){
    perc.cont$Chr[i] = paste(strsplit(perc.cont$Chr,"_")[[i]][1],(strsplit(perc.cont$Chr,"chr")[[i]][2]),sep="_")
}
write.csv(perc.cont,f.cont,row.names=F)
