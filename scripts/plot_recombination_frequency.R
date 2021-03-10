#!/usr/bin/env Rscript

##############################################################
#  script: plot_recombination_frequency.R
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.10.12
#  description: plot recombination frequency based on event frequency sliding window  file
#  example: Rscript --vanilla plot_recombination_frequency.R -i pooled.events.recombination_frequency.W10000S2000.txt -c ref.centromere.gff --output output.pdf
##############################################################

library("ggplot2")
library("optparse")

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "input genotype file name", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "out.pdf", 
              help = "output plot file name [default= %default]", metavar = "character"),
  make_option(c("-c", "--centromere"), type = "character", default = "NULL", 
              help = "the centromere annotation in gff format", metavar = "character")
); 

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("the input file is missing", call. = FALSE)
}

data <- read.table(opt$input, sep = "\t", header = FALSE)
colnames(data) <- c("chr", "pos", "event_count", "event_type")

centromere_data <- read.table(opt$centromere, sep = "\t", header = FALSE)
colnames(centromere_data) <- c("chr", "tag", "feature", "start", "end", "score",
                                "strand", "phase", "attributes")

chr_list <- as.vector(unique(data$chr))
chr_list <- head(chr_list, -1)
ymax <- max(data$event_count, na.rm = TRUE)

color_palette <- c("CO" = "darkorange", "GC" = "mediumpurple1", "NCO" = "dodgerblue")

pdf(file = opt$output, height = 4, width = 12)

for (c in chr_list) {
  subdata <- subset(data, chr == c)
  centromere_subdata <- subset(centromere_data, chr == c)
  print(ggplot(subdata, aes(x = pos, y = event_count, color = event_type)) + 
    geom_line(size = 0.5) + 
    scale_color_manual(values = color_palette) +
    geom_rect(data = centromere_subdata, aes(xmin = start, xmax = end, ymin = 0, ymax = ymax),
              color = "red",
              fill = "red",
              alpha = 1,
              inherit.aes = FALSE) +
    scale_x_continuous(name=paste(c, " (kb)"), 
        breaks = seq(0, tail(subdata$pos,1), 100000), 
 	labels = seq(0, tail(subdata$pos,1)/1000, 100),
	minor_breaks = seq(0, tail(subdata$pos, 1), 50000)) +
    scale_y_continuous(name = "Event count") + 
    # ggtitle(paste0(c)) +
    facet_wrap(~event_type, ncol = 1, scale = "free_x") +
    theme_bw())
  #ggsave(output, height = 4, width = 1.0*tail(subdata$pos,1)/50000, limitsize = FALSE)
}




