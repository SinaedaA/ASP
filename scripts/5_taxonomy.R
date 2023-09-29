#!/usr/bin/env Rscript

#### ENV SET-UP AND ARGUMENTS ####
##### Loading packages silently #####
# Check if one of .libPaths() is writable for the user, and if not, create a writable package library in $HOME (~)
source("~/Rfunctions/check_libPaths.R")
writable_path <- check_libPaths(verbose = FALSE)
# Install/load pacman.
suppressPackageStartupMessages(if (!require(pacman)) {
    install.packages("pacman", lib = writable_path)
})
needed_packs <- c("BiocManager", "tidyverse", "plyr", "ggplot2", "Biostrings", "cli", "argparser", "tibble", "this.path", "DECIPHER") # , "here"
github_packs <- c("benjjneb/dada2")
to_install <- needed_packs[!needed_packs %in% pacman::p_library()]
message("Installing and/or loading required packages...")
if (length(to_install) != 0) {
    pacman::p_install(to_install, character.only = TRUE, try.bioconductor = TRUE, path = writable_path)
}
pacman::p_load(needed_packs, character.only = TRUE)
pacman::p_load_gh(github_packs)

script_dir <- this.dir()
working_dir <- getinitwd()

##### Parse arguments with arg_parser #####
cli_h1("Parsing input")
parser <- arg_parser("Taxonomy assignment of denoised amplicon sequencing data")
parser <- add_argument(parser, "sequence_table", help = "Absolute path to directory where representative sequence count tables are located (seqtab_nochim.tsv)")
parser <- add_argument(parser, "output_directory", help = "Absolute path to directory where the output should be written.")
parser <- add_argument(parser, "taxonomy", help = "Absolute path to the appropriate taxonomy training set (downloaded from DECIPHER website).")
argv <- parse_args(parser)

#### Load data ####
count_table <- read.table(argv$sequence_table)
load(argv$taxonomy)

#### Load DNA from count_table colnames ####
dna <- DNAStringSet(getSequences(colnames(count_table))) # Create a DNAStringSet from the ASVs

#### Match to taxonomy training set ####
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors

#### Get taxonomy data and reformat to look like output of assignTaxonomy (from dada2) ####
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(colnames(count_table))

#### Make taxonomy table ####
taxonomy <- as.data.frame(taxid) %>% 
    dplyr::mutate(ASV = paste0("ASV", seq(nrow(taxid)))) %>%
    rownames_to_column("Seq") %>%
    column_to_rownames("ASV") %>%
    relocate("domain", "phylum", "class",  "order",  "family", "genus",  "species", "Seq")
## make correspondance table between sequences and ASV numbers
asv2seq <- setNames(as.data.frame(cbind(rownames(taxonomy), taxonomy$Seq)), nm = c("ASV", "Seq"))
## change the column names of count_table
count_table <- count_table %>%
  rename_with(~coalesce(asv2seq$ASV[match(., asv2seq$Seq)], .))

write.table(count_table, file = paste0(argv$output_directory, "/count_table.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(asv2seq, file = paste0(argv$output_directory, "/asv2seq.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(taxonomy, file = paste0(argv$output_directory, "/taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)