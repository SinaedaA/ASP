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
needed_packs <- c("BiocManager", "tidyverse", "plyr", "ggplot2", "Biostrings", "cli", "argparser", "tibble", "this.path", "DECIPHER", "phylotools") # , "here"
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

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

##### Parse arguments with arg_parser #####
cli_h1("Parsing input")
parser <- arg_parser("Taxonomy assignment of denoised amplicon sequencing data")
parser <- add_argument(parser, "sequence_table", help = "Absolute path to directory where representative sequence count tables are located (seqtab_nochim.tsv)")
parser <- add_argument(parser, "output_directory", help = "Absolute path to directory where the output should be written.")
parser <- add_argument(parser, "taxonomy", help = "Absolute path to the appropriate taxonomy training set (downloaded from DECIPHER website).")
argv <- parse_args(parser)

#### Load data ####
cli_h1("Loading count table")
count_table <- read.table(argv$sequence_table)
print(length(colnames(count_table)))
load(argv$taxonomy)

#### Make output directory ####
cli_h1("Checking existence of output directory...")
working_dir <- getinitwd()
outpath <- paste0(working_dir, "/", argv$output_directory)
if (!dir.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
    cli_alert_info("Output directory {outpath} doesn't exist, creating it.")
} else {
    cli_alert_info("Output directory {outpath} already exists.")
}

#### Load DNA from count_table colnames ####
cli_h1("Extracting DNA sequences")
dna <- DNAStringSet(getSequences(colnames(count_table))) # Create a DNAStringSet from the ASVs

#### Match to taxonomy training set ####
cli_h1("Matching sequences to training set")
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=TRUE) # use all processors

#### Get taxonomy data and reformat to look like output of assignTaxonomy (from dada2) ####
cli_h1("Formatting output tables.")
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

cli_h1("Writing output tables...")
write.table(count_table, file = paste0(outpath, "/count_table.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
cli_alert_success("Count table written to {outpath}/count_table.tsv : ASV x Sample")
write.table(asv2seq, file = paste0(outpath, "/asv2seq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cli_alert_success("ASV to SEQ correspondance written to {outpath}/asv2seq.tsv : ASV x SEQ")
write.table(taxonomy, file = paste0(outpath, "/taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
cli_alert_success("Taxonomy table written to {outpath}/taxonomy.tsv : ASV -> domain -> phylum -> class -> order -> family -> genus -> species -> Seq")
writeFasta(setNames(as_tibble(asv2seq), nm = c("name", "seq")), paste0(outpath, "ASV.fasta"))
cli_alert_success("ASV fasta file written to {outpath}/ASV.fasta : >ASVxxx \n sequence")