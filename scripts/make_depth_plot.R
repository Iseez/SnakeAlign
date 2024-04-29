library(data.table)
library(zoo)
library(ggplot2)

if( getDTthreads() < 10 || getDTthreads() > 10){
  setDTthreads(threads = 10)
}else{print(paste("WARNING: data.table package is running with ", getDTthreads(), " threads.", sep=''))}

plus_10 <- function(y){
    n = length(y)
    pt = 0
    for(x in y){
        if(x > 10){
            pt = pt + 1
        }
    }
    return(format(round((pt/n)*100, 4),nsmall=4))
}

FILE = snakemake@input[[1]]

SIZE = snakemake@config[["size"]]
STEP = snakemake@config[["step"]]

PerContig = snakemake@output[[1]]
Windows = snakemake@output[[2]]
PerRefPATH = snakemake@output[[3]]
imgdpth = snakemake@output[[4]]
imgdpth.m = snakemake@output[[5]]

depth <- fread(FILE)
names(depth) <- c("Chr", "locus", "depth")

depth[,`:=`(Ref = vapply(strsplit(Chr,"_"), `[`, 1, FUN.VALUE=character(1)) )]
depth2plot <- depth[, .(window.start = rollapply(locus, width=SIZE, by=STEP, FUN=min, align="left"),window.end = rollapply(locus, width=SIZE, by=STEP, FUN=max, align="left"),coverage.mean = rollapply(depth, width=SIZE, by=STEP, FUN=mean, align="left")), .(Chr)]
write.csv(depth2plot,Windows,row.names=FALSE)

depth2plot[,`:=`(Ref = vapply(strsplit(Chr,"_"), `[`, 1, FUN.VALUE=character(1)) )]
depth2plot[,`:=`(chr = vapply(strsplit(Chr,"_"), `[`, 5, FUN.VALUE=character(1)) )]
depth2plot[,`:=`(Chr = paste(Ref,chr,sep="_chr") )]

nomt <- depth2plot[chr!="mt" & chr!="2m"]
mt <- depth[grep("mt",Chr)]

upr <- ceiling(quantile(nomt$coverage.mean,probs=0.999))
uprm <- ceiling(max(mt$depth))

p <- ggplot(nomt, aes(x=window.end, y=coverage.mean, colour=Ref)) +
    geom_point(size = 1,alpha=0.5,shape=16) +
    scale_x_continuous(name="Genomic Position in bp (winSize=1Kb)") +
    scale_y_continuous(name="Average Coverage Depth", limits=c(0, upr), breaks = seq(0,upr,25)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_grid(. ~ Chr, space="free_x", scales = "free_x") +
    theme(panel.spacing.x = grid::unit(0, "cm"),panel.border = element_rect(colour = "grey",size = 0.1),panel.ontop = FALSE) +
    theme(text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x=element_blank(),
    strip.background = element_blank(), strip.text.x = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 7)))

ggsave(imgdpth,plot = p, width = 20, height = 12, units = "cm",bg = "white", dpi = 300)

p <- ggplot(mt, aes(x=locus, y=depth, colour=Ref)) +
    geom_point(shape = 16, size = 1,alpha=0.5) +
    scale_x_continuous(name="Genomic Position in bp (winSize=1bp)") +
    scale_y_continuous(name="Average Coverage Depth", limits=c(0, uprm), breaks = seq(0,uprm,200)) +
    scale_shape_discrete(name  ="Reference") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_grid(. ~ Chr, space="free_x", scales = "free_x") +
    theme(panel.spacing.x = grid::unit(0, "cm"),panel.border = element_rect(colour = "grey",size = 0.05)) +
    theme(text = element_text(size=10), axis.text.x = element_blank(),axis.ticks.x=element_blank(),
    strip.background = element_blank(), strip.text.x = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 7)))

ggsave(imgdpth.m,plot = p, width = 20, height = 12, units = "cm",bg = "white", dpi = 300)

gc()

depth_contig <- depth[,.(Mean = mean(as.double(depth)),Min = min(as.double(depth)),
  Max = max(as.double(depth)),l_quantile = quantile(as.double(depth),probs=0.25),Median = median(as.double(depth)),
  u_quantile = quantile(as.double(depth),probs=0.75),Above_10X = plus_10(as.double(depth))),.(Chr)]
for(i in 1:length(depth_contig$Chr)){
    depth_contig$Chr[i] = paste(strsplit(depth_contig$Chr,"_")[[i]][1],(strsplit(depth_contig$Chr,"_")[[i]][5]),sep="_")
}
write.csv(depth_contig,PerContig,row.names=FALSE)

genome_Stats <- depth[,.(Reads = sum(as.double(depth)),Mean = mean(as.double(depth)),Min = min(as.double(depth)),
  Max = max(as.double(depth)),l_quantile = quantile(as.double(depth),probs=0.25),Median = median(as.double(depth)),
  u_quantile = quantile(as.double(depth),probs=0.75),Above_10X = plus_10(as.double(depth))),.(Ref)]
write.csv(genome_Stats,PerRefPATH,row.names=FALSE)
