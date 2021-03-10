#!/usr/bin/env Rscript
##############################################################
#  script: plot_parental_allele_frequency_in_tetrads.R
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.09.22
#  description: plot genome-wide parental allele frequency in tetrads
#  example: Rscript --vanilla --slave plot_parental_allele_frequency_in_tetrads.R --input input.parental_allele_frequency.txt.gz --color color_scheme.txt --output output.parental_allele_frequency.plot.pdf  
##############################################################

library("optparse")
library("ggplot2")

option_list <- list(
  make_option(c("--input"), type = "character", default = NULL, 
              help = "input genotype file name", metavar = "character"),
  make_option(c("--color_scheme"), type = "character", default = NULL, 
                help = "the color scheme used for the plot", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL, 
              help = "output file name", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

data <- read.table(opt$input, header = TRUE, sep = "\t", na.strings = "")

# read color scheme
color_scheme <- read.table(opt$color_scheme, header = FALSE, sep = "\t", comment.char = "", na.strings = "")
colnames(color_scheme) <- c("genotype", "color")
color_palette <- as.character(color_scheme$color)
names(color_palette) <- as.character(color_scheme$genotype)

chr_list <- as.vector(unique(data$chr))
pdf(file = opt$output, height = 3, width = 8)

for (c in chr_list) {
    subdata <- subset(data, chr == c)
    print(ggplot(subdata, aes(x = pos, y = freq, color = genotype, alpha = 0.5)) + 
    	geom_point(shape = 3, size = 0.5) +
    	scale_color_manual(values = color_palette) +
	scale_x_continuous(name=paste(c, " (kb)"), breaks=seq(0,tail(subdata$pos,1),50000), labels=seq(0,tail(subdata$pos,1)/1000,50)) + 
	scale_y_continuous("Parental allele frequency") +
	facet_wrap(~chr, ncol = 1, scale = "free_x") +
    	theme_bw() +
	theme(axis.text.x = element_text(angle=45,hjust=1)))
}

dev.off()


