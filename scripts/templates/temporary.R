# count table
freqtable <- read.table("~/Documents/bioinfo/decrypt/ASP/fungi/count_table.tsv", sep = "\t", header = TRUE, comment.char = "")
colnames(freqtable) <- gsub("\\.", "-", colnames(freqtable))
freqtable$taxonomy <- gsub(".__", "", freqtable$taxonomy)
# metadata
columns <- c("Sample_Name", "Short_Sample", "Sample_ID", "Organism", "Genotype", "Compartment", "Timepoint", "Abbrev", "Summary")
metadata <- setNames(read.table("~/Documents/bioinfo/decrypt/ASP/fungi/sample-metadata-withEmpty.tsv", header = TRUE, sep = "\t", comment.char = ""), nm = columns)
if (!identical(sort(colnames(freqtable)[-1]), sort(metadata$Sample_Name))) {
    cli_alert_danger("The colnames of the count table (sample names) do not correspond to the Sample_Name column in metadata")
    stop()
}

freqtable <- freqtable %>%
    separate(col = taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
    replace(is.na(.), "unassigned") %>%
    mutate_all(list(~ str_remove_all(., "D_")))

metadata <- subset(metadata, metadata$Abbrev == "Hv")
metadata$Genotype[!metadata$Genotype == "BI1"] <- "WT"
freqtable <- freqtable %>%
    dplyr::select(1:7, metadata$Sample_Name)

taxon_data_list <- MakeTaxonDataList(freqtable, taxons)
filtered_list <- FilterAbs(taxon_data_list, metadata, taxons, Nreads = 2000)
transformed_list <- Transform(filtered_list, taxons, Filter = "TotalAbundance")
alpha_list <- GetAlphaList(metadata, filtered_list, taxons, index = "shannon")

## Boxplots for relative abundance differences:
transformed <- transformed_list[["genus"]]
taxa_transf <- transformed %>%
    tibble::rownames_to_column("Sample_Name") %>%
    dplyr::left_join(metadata[c("Sample_Name", "Compartment", "Genotype", "Timepoint")], by = "Sample_Name") %>%
    tidyr::pivot_longer(-c("Sample_Name", "Compartment", "Genotype", "Timepoint"), names_to = "Taxon", values_to = "Abundance") %>%
    dplyr::mutate(Compartment = factor(Compartment, levels = c("Roots", "Rhizoplane", "Rhizosphere", "Soil")))
kw_test <- taxa_transf %>%
    group_by(Compartment, Timepoint, Taxon) %>%
    group_map(~ KW_wrapper(.y, .x$Abundance, .x$Genotype)) %>%
    bind_rows()
taxa_transf <- taxa_transf %>%
    dplyr::left_join(kw_test, by = c("Compartment", "Timepoint", "Taxon"))
# skipping the looping over compartments, and only considering roots
Comp <- "Roots"
Comp <- "Rhizosphere"
data <- dplyr::filter(taxa_transf, Compartment == Comp)
tax_data <- data %>%
    dplyr::filter(p_value <= 0.05) %>%
    dplyr::filter(Taxon != "Fusarium" & Taxon != "Ascobolus")
up_tax_data <- tax_data %>%
    dplyr::group_by(Genotype, Taxon) %>%
    mutate(MeanAbundance = mean(Abundance)) %>% 
    dplyr::filter(mean())
taxa <- tax_data %>%
    dplyr::select(Taxon) %>%
    unique()
Colors = myColors

## would need to order them based on whether they are up or down. Maybe separate them, then patch
## add fusarium separately
## add label for timepoint at which this happens

p <- tax_data %>%
    ggplot(aes(y = factor(Genotype, levels = c("BI1", "WT")), x = Abundance, color = Compartment)) +
    geom_boxplot(outlier.shape = NA, width = .5) + # outlier.size = 3, outlier.colour = "black", ) ++
    coord_cartesian(xlim = quantile(tax_data$Abundance, c(0, 0.99)))+
    geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
    facet_grid(Taxon ~ Compartment, scales = "free") +
    scale_color_manual(values = Colors) +
    theme_bw() +
    theme(
        strip.text = element_text(size = 20), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)
    ) +
    guides(x = guide_axis(angle = 45)) +
    labs(y = "Genotype", x = "Relative abundance") 