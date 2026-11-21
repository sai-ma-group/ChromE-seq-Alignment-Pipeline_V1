#!/usr/bin/Rscript
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
args <- commandArgs()
#print(args)

dir <- args[6]
Name <- args[7]
out <- args[8]

# read data
Counts <- data.table::fread(paste(dir, Name, sep="/"), sep="\t")
df <- Counts %>% group_by(V5) %>% count() %>%  arrange(desc(n))


pdf(out)
plot(log10(df$n), xlab="Barcode rank", ylab = "log10 (contacts)", main = "No. contacts per cell", col="darkblue", pch=16)
garbage <- dev.off()
