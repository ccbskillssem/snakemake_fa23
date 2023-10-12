# This script is  derived from the BAGEP pipeline, where I've replaced the heatmap generation function with a simplified version for this workshop.
# https://github.com/idolawoye/BAGEP
# Many thanks to the BAGEP team for their easy-to-parse code!
# This script supports the vcf heatmap generation rule that I borrowed to demonstrate how R scripts can be incorporated into snakemake.

library('vcfR')
library('argparse')

parser <- ArgumentParser(description='VCF visualisation script')
parser$add_argument("vcf", type='character', action="store", help="Input VCF file")
parser$print_help()
args <- commandArgs(trailingOnly = TRUE)

if (length(args)>2) {
  stop("Only two inputs must be supplied: VCF file and plot outfile path", call.=FALSE)
} else if (length(args)==3) {
  noquote(". . . Checking input files")
}
vcf_file <- args[[1]]
vcf <- read.vcfR(vcf_file, verbose=FALSE, nrows = 500)

genotype <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
rownames(genotype) <- 1:nrow(genotype)
genotype[!rowSums(!is.finite(genotype)),]
genotype[!is.finite(genotype)] <- 0

# this creates a ggplot object, but ggsave doesn't reliably play well with RScript execution
# so we save using the pdf creation + print option
pdf(args[[2]])
vcf_heatmap <- heatmap.bp(genotype)
print(vcf_heatmap)
dev.off()