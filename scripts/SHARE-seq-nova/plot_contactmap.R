#!/usr/bin/Rscript
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(data.table))
args <- commandArgs()
options(warn = -1)
#print(args)

# Rscript $myPATH/plot_contactmap.R $dir/fastqs/ $Name.$Species.contacts.final.gz $Species 500000 $Name.$Species --save
# /mnt/Scripts/SHARE-seq-nova/plot_contactmap.R ./ DNA-scHiC-SpeciesMix-1K-2024-05-22.hg38.contacts.final.gz hg38 500000 DNA-scHiC-SpeciesMix-1K-2024-05-22.hg38 --save

dir <- args[6]
Name <- args[7]
Species <- args[8]
resolution <- as.numeric(args[9])
out <- args[10]

# read data

splitAndFetch <- function(vec, delim, part){
if (length(part) == 1) {
        sapply(strsplit(as.character(vec), delim, fixed = TRUE), 
            "[[", part)
    }
    else {
        sapply(strsplit(as.character(vec), delim, fixed = TRUE), 
            function(x) paste(x[part], collapse = delim))
    }
}
contacts2MX <-function(file, genome, resolution){
	    df <- data.table::fread(file, sep="\t")
	    BCs <- unique(df$V5)
	    df.sub <- df %>% filter(V5 %in% BCs)
	    if (genome == "mm10"){
	      require(BSgenome.Mmusculus.UCSC.mm10)
	      chrs <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)[1:22]
	      chrom.size <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)[chrs]
	    } else if (genome == "hg38"){
	      require(BSgenome.Hsapiens.UCSC.hg38)
	      chrs <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)[1:25]
	      chrom.size <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[chrs]
	    } else if (genome == "hg19"){
	    require(BSgenome.Hsapiens.UCSC.hg19)
	    chrs <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)[1:25]
	    chrom.size <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chrs]
	    }

bins <- tileGenome(chrom.size, tilewidth=resolution, cut.last.tile.in.chrom=T)

gr1 <- makeGRangesFromDataFrame(df.sub, seqnames.field = "V1", start.field = "V2", end.field = "V2")
gr2 <- makeGRangesFromDataFrame(df.sub, seqnames.field = "V3", start.field = "V4", end.field = "V4")

ov1 <- findOverlapPairs(gr1, bins)
df.sub$bin1 <- paste(seqnames(ov1@second), end(ov1@second), sep=":")

ov2 <- findOverlapPairs(gr2, bins)
df.sub$bin2 <- paste(seqnames(ov2@second), end(ov2@second), sep=":")

summarized <- df.sub %>% group_by(V5, bin1, bin2) %>% count() 
BCs.identified <- unique(summarized$V5)
bcmap <- data.frame(BCs=BCs.identified)
bcmap$cell_id <- 1:nrow(bcmap)

contact.mx <- data.frame(cell_id = match(summarized$V5,bcmap$BCs), 
                    chrom1 = splitAndFetch(vec=summarized$bin1,delim=":",part=1),
                    pos1 = splitAndFetch(vec=summarized$bin1,delim=":",part=2),
                    chrom2 = splitAndFetch(vec=summarized$bin2,delim=":",part=1),
                    pos2 = splitAndFetch(vec=summarized$bin2,delim=":",part=2),
                    count = summarized$n)
return(list(bcmap, contact.mx, chrom.size))
}

input <- paste(dir, "/", Name, sep="")
# input
contact.mx <- contacts2MX(input, Species, resolution)
# head(contact.mx)

bcmap <- contact.mx[[1]]
onecell <- contact.mx[[2]]
chrom.size <- contact.mx[[3]]

gr1 <- makeGRangesFromDataFrame(onecell, seqnames.field = "chrom1", start.field = "pos1", end.field = "pos1")
gr2 <- makeGRangesFromDataFrame(onecell, seqnames.field = "chrom2", start.field = "pos2", end.field = "pos2")
bins <- tileGenome(chrom.size, tilewidth=resolution, cut.last.tile.in.chrom=T)

bins$rank <- 1:length(bins)
ov1 <- findOverlaps(gr1, bins)
onecell$bin1 <- subjectHits(ov1)
ov2 <- findOverlaps(gr2, bins)
onecell$bin2 <- subjectHits(ov2)

temp <- onecell
temp$bin1 <- onecell$bin2
temp$bin2 <- onecell$bin1
df <- rbind(onecell, temp)
#df <- df[df$chrom1 %in% chrs & df$chrom2 %in% chrs, ]

png(paste0(dir, "/", out, ".contactmap.png"), width=4.25, height=3, units = "in", res=1200)
ggplot(df, aes(x=bin1, y=bin2, fill=count)) +
  geom_tile() +
  theme_classic() +
  ggtitle(Species) +
  scale_y_reverse() +
 # xlab(paste(chrs[], collapse = '')) +
 # ylab(paste(chrs[], collapse = ''))
  xlab("chrs") +
  ylab("chrs") 
  #scale_fill_viridis(option = "C")
garbage <- dev.off()

## observed contact map
#Mx <- reshape2::dcast(df[, 6:8], bin1 ~ bin2)
Mx <- dcast(df[, 6:8], bin1 ~ bin2, value.var = "bin2", fun.aggregate = sum)
rownames(Mx) <- Mx$bin1
Mx <- Mx[, 2:ncol(Mx)]
Mx <- as.matrix(Mx)

## expected contact map
temp <- proxy:::dist(1:nrow(Mx), 1:nrow(Mx))

## correlation
corr <- cor(Mx*temp/100, Mx*temp/100, method = "spearman")
corr[is.na(corr)] <-0
limits <- quantile(corr, prob=c(0.0, 0.95))
mat_breaks <- seq(limits[1], limits[2], length.out = 50)

library(pheatmap)
colors <- c("#053061", "#2971B1", "#6AACD0", "#C1DDEB", "#F7F7F7", "#FACDB5", "#E58267", "#BB2933", "#67001F")
png(paste0(dir, "/", out, ".corrmap.png"), width=4, height=3, units = "in", res=1200)
pheatmap(corr, cluster_col = F, cluster_rows = F, show_colnames = F, show_rownames = F,
         scale="none",
         main = "Auto Correlation", 
         # color = rev(viridis(length(mat_breaks) - 1, option = "A")), # inferno heat.colors  viridis
         # color = rev(heat.colors(length(mat_breaks) - 1)), # inferno heat.colors  viridis
         # color = coolerColors(),
         color = colorRampPalette(colors)(length(mat_breaks)), # wolfgang_extra brewer_celsius brewer_yes
         breaks = mat_breaks)
garbage <- dev.off()
