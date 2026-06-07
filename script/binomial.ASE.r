#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript binomial.ASE.r <input.tested.sites> <output.ASE.sites>")
}

infile <- args[1]
outfile <- args[2]

df <- fread(infile, sep = "\t", header = TRUE)

df[, p := mapply(function(ref, alt) {
  binom.test(ref, ref + alt, p = 0.5)$p.value
}, refCount, altCount)]

df[, FDR := p.adjust(p, method = "BH")]
df[, effect := abs(0.5 - ref_ratio)]

fwrite(df, outfile, sep = "\t", quote = FALSE)
