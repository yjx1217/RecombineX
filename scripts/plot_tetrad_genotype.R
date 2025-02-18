#!/usr/bin/env Rscript

##############################################################
#  script: plot_tetrad_genotype.R
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2025.02.18
#  description: plot tetrad genotype based on the markers
#  example: Rscript --vanilla --slave plot_tetrad_genotype.R --input input.genotype.lite.raw.for_genotype_plotting.txt.gz --color color_scheme.txt --genome1_tag genome1_tag --genome2_tag genome2_tag --coordinate_genome_fai genome1_tag.genome.raw.fa.fai --coordinate_genome_centromere_gff genome1_tag.centromere.gff --output output.plot.pdf (--query_chr query_chr --query_start query_start --query_end query_end --query_flanking 1000) 
##############################################################

library("optparse")
library("ggplot2")

option_list <- list(
  make_option(c("--input"), type = "character", default = NULL, 
              help = "input genotype file name", metavar = "character"),
  make_option(c("--genome1_tag"), type = "character", default = NULL, 
              help = "the genome tag for genome1", metavar = "character"),
  make_option(c("--genome2_tag"), type = "character", default = NULL, 
              help = "the genome tag for genome2", metavar = "character"),
  make_option(c("--coordinate_genome_fai"), type = "character", default = NULL, 
              help = "the .fai file for the genome used to for set up genome coordinates", metavar = "character"),
  make_option(c("--coordinate_genome_centromere_gff"), type = "character", default = NULL, 
              help = "the centromere GFF3 file for the genome used to for set up genome coordinates", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL, 
              help = "output plot file name", metavar = "character"),
  make_option(c("--color_scheme"), type = "character", default = NULL, 
              help = "the color scheme used for the plot", metavar = "character"),
  make_option(c("--query_chr"), type = "character", default = NULL,
              help = "query chromosome", metavar = "character"),
  make_option(c("--query_start"), type = "integer", default = NULL,
              help = "the start coordinate of the query region", metavar = "character"),
  make_option(c("--query_end"), type = "integer", default = NULL,
              help = "the end coordinate of the query region", metavar = "character"),
  make_option(c("--query_flanking"), type = "integer", default = 0,
              help = "the flanking radius of the query region", metavar = "character")
); 

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("the genotype input file is missing", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("the output plot file is missing", call. = FALSE)
}
if (is.null(opt$color_scheme)) {
  print_help(opt_parser)
  stop("the color scheme file is missing", call. = FALSE)
}

# read genotype

genome1_tag <- opt$genome1_tag
genome2_tag <- opt$genome2_tag

genotype <- read.table(opt$input, header = TRUE, sep = "\t", na.strings = "")

genotype$gamete_id <- genotype$spore_id
genotype$gamete_genotype <- factor(genotype$spore_genotype)


# read color scheme
color_scheme <- read.table(opt$color_scheme, header = FALSE, sep = "\t", na.strings = "")
colnames(color_scheme) <- c("genome_tag", "color")
color_palette <- as.character(color_scheme$color)
names(color_palette) <- as.character(color_scheme$genome_tag)


# setup the genome and chromosomes
genome_size <- read.table(opt$coordinate_genome_fai, header = FALSE, sep = "\t")
colnames(genome_size) <- c("chr", "genome_size", "offset", "linebases", "linewidth")
genome_size_dict <- genome_size$genome_size
names(genome_size_dict) <- genome_size$chr

if (!is.null(opt$coordinate_genome_centromere_gff)){
   centromere <- read.table(opt$coordinate_genome_centromere_gff, header = FALSE, sep = "\t")
   colnames(centromere) <- c("chr", "ref", "type", "start", "end", "score", "strand", "phase", "info")
}

if (!is.null(opt$query_chr)){
   if (is.null(opt$query_start) || is.null(opt$query_end)) {
     stop("query_chr, query_start, and query_end need to be specified in the same time", call. = FALSE)
   } else if (opt$query_start > opt$query_end) {
     stop("Error: query_start > query_end", call. = FALSE)
   } else if (is.null(opt$query_flanking)) {
     opt$query_flanking <- 0;
   }
   plot_chr <- opt$query_chr
   plot_start <- opt$query_start - opt$query_flanking 
   plot_end <- opt$query_end + opt$query_flanking 
   
   genotype <- subset(genotype, chr == opt$query_chr)
   
   if (opt$query_start - opt$query_flanking <= 0) {
     plot_start <- 1
   }
   if (opt$query_end + opt$query_flanking > genome_size_dict[plot_chr]) {
     plot_end <- genome_size_dict[plot_chr]
   }
}


# re-order the chromosomes

chr_list <- factor(genome_size$chr, levels = unique(genome_size$chr))
genome_size$chr <-factor(chr_list, levels = rev(chr_list))

gamete_indice <- c("a", "b", "c", "d")

chr_gamete_list <- list()
for (c in genome_size$chr) { 
  chr_gamete_list <- append(chr_gamete_list, paste0(c, ".", gamete_indice))
}

genome_size.expanded <- genome_size[rep(row.names(genome_size), each = 4), 1:2]
genome_size.expanded$gamete_id <- rep(gamete_indice, length(genome_size$chr))
genome_size.expanded$chr_gamete <- paste0(genome_size.expanded$chr, ".", genome_size.expanded$gamete_id)
genome_size.expanded <- genome_size.expanded[match(chr_gamete_list, genome_size.expanded$chr_gamete), ]
genome_size.expanded$chr_gamete <- factor(chr_gamete_list, levels = rev(chr_gamete_list))

genotype$chr_gamete <- paste0(genotype$chr, ".", genotype$gamete_id)

if (!is.null(opt$coordinate_genome_centromere_gff)){
   centromere.expanded <- centromere[rep(row.names(centromere), each = 4), 1:9]
   centromere.expanded$gamete_id <- rep(gamete_indice, length(centromere$chr))
   centromere.expanded$chr_gamete <- paste0(centromere.expanded$chr, ".", centromere.expanded$gamete_id)
}

if (!is.null(opt$query_chr)){
   genome_size.expanded <- subset(genome_size.expanded, chr == opt$query_chr)
   if (!is.null(opt$coordinate_genome_centromere_gff)) {
       centromere.expanded <- subset(centromere.expanded, chr == opt$query_chr)
   }
}

if (!is.null(opt$coordinate_genome_centromere_gff) && is.null(opt$query_chr)) {
   ggplot(data = genome_size.expanded, aes(x = chr_gamete, y = genome_size)) +
       geom_bar(stat = "identity", width = 0.5, fill = "grey90", color = "black") +
       geom_segment(data = genotype, aes(x = chr_gamete, xend = chr_gamete, 
       y = adjusted_marker_start, yend = adjusted_marker_end, 
       color = gamete_genotype), linewidth = 0.8) +
       scale_color_manual(values = color_palette) +
       geom_point(data = centromere.expanded, aes(x = chr_gamete, y = (start + end) / 2), 
       shape = 21, size = 1.5, fill = "white") + 
       scale_y_continuous(name = "Size (bp)") +
       scale_x_discrete(name = "Chromosome") +
       coord_flip() + 
       theme_classic()
       ggsave(filename = opt$output, device = "pdf", width = 10, height = 8, limitsize = FALSE)
} else if (is.null(opt$coordinate_genome_centromere_gff) && is.null(opt$query_chr)) {
   ggplot(data = genome_size.expanded, aes(x = chr_gamete, y = genome_size)) +
       geom_bar(stat = "identity", width = 0.5, fill = "grey90", color = "black") +
       geom_segment(data = genotype, aes(x = chr_gamete, xend = chr_gamete, 
       y = adjusted_marker_start, yend = adjusted_marker_end, 
       color = gamete_genotype), linewidth = 0.8) +
       scale_color_manual(values = color_palette) + 
       scale_y_continuous(name = "Size (bp)") +
       scale_x_discrete(name = "Chromosome") +
       coord_flip() + 
       theme_classic()
       ggsave(filename = opt$output, device = "pdf", width = 10, height = 8, limitsize = FALSE)
} else if (!is.null(opt$coordinate_genome_centromere_gff) && !is.null(opt$query_chr)) {
   ggplot(data = genome_size.expanded, aes(x = chr_gamete, y = genome_size)) +
       geom_bar(stat = "identity", width = 0.5, fill = "grey90", color = "black") +
       geom_segment(data = genotype, aes(x = chr_gamete, xend = chr_gamete, 
       y = raw_marker_start - 2, yend = raw_marker_end + 2, 
       color = gamete_genotype), linewidth = 5) +
       scale_color_manual(values = color_palette) + 
       geom_point(data = centromere.expanded, aes(x = chr_gamete, y = (start + end) / 2), 
       shape = 21, size = 1.5, fill = "white") + 
       scale_y_continuous(name = "Size (bp)") +
       scale_x_discrete(name = "Chromosome") +
       coord_flip(ylim = c(plot_start, plot_end)) + 
       theme_classic()
       ggsave(filename = opt$output, device = "pdf", width = 6, height = 4, limitsize = FALSE)
} else {
   ggplot(data = genome_size.expanded, aes(x = chr_gamete, y = genome_size)) +
       geom_bar(stat = "identity", width = 0.5, fill = "grey90", color = "black") +
       geom_segment(data = genotype, aes(x = chr_gamete, xend = chr_gamete, 
       y = raw_marker_start - 2, yend = raw_marker_end + 2, 
       color = gamete_genotype), linewidth = 5) +
       scale_color_manual(values = color_palette) + 
       scale_y_continuous(name = "Size (bp)") +
       scale_x_discrete(name = "Chromosome") +
       coord_flip(ylim = c(plot_start, plot_end)) + 
       theme_classic()
       ggsave(filename = opt$output, device = "pdf", width = 6, height = 4, limitsize = FALSE)
}
