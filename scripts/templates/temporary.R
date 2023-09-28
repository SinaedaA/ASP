source("~/Documents/projects/ASP/scripts/data_analysis_options.Rdata")
# count table
#freqtable <- read.table("~/Documents/bioinfo/decrypt/ASP/fungi/count_table.tsv", sep = "\t", header = TRUE, comment.char = "")
freqtable <- read.table("~/Documents/projects/ASP/fungi/count_table.tsv", sep = "\t", header = TRUE, comment.char = "")
colnames(freqtable) <- gsub("\\.", "-", colnames(freqtable))
freqtable$taxonomy <- gsub(".__", "", freqtable$taxonomy)
# metadata
columns <- c("Sample_Name", "Short_Sample", "Sample_ID", "Organism", "Genotype", "Compartment", "Timepoint", "Abbrev", "Summary")
metadata <- setNames(read.table("~/Documents/projects/ASP/fungi/sample-metadata-withEmpty.tsv", header = TRUE, sep = "\t", comment.char = ""), nm = columns)
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
outdir <- "~/Documents/projects/ASP/fungi/data_analysis/"

taxon_data_list <- MakeTaxonDataList(freqtable, taxons)
filtered_list <- FilterAbs(taxon_data_list, metadata, taxons, Nreads = 2000, Outdir = outdir)
transformed_list <- Transform(filtered_list, taxons, Filter = "TotalAbundance")
alpha_list <- GetAlphaList(metadata, filtered_list, taxons, index = "shannon")
alpha_plots <- AlphaPlotList(alpha_list, taxons, myColors, index = "shannon")

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
# Comp <- "Rhizosphere"
data <- dplyr::filter(taxa_transf, Compartment == Comp)
tax_data <- data %>%
    dplyr::filter(p_value <= 0.05) %>%
    dplyr::filter(Taxon != "Fusarium")
fusarium_data <- data  %>% 
    dplyr::filter(p_value <= 0.05 & Taxon == "Fusarium")
low_tax_data <- dplyr::filter(tax_data, Taxon %in% c("Ascobolus", "Cladorrhinum", "Iodophanus", "Slopeiomyces", "Trichocladium"))
mid_tax_data <- dplyr::filter(tax_data, Taxon %in% c("Bipolaris", "Dactylella", "Penicillium", "Periconia", "Schizothecium"))
taxa <- tax_data %>%
    dplyr::select(Taxon) %>%
    unique()
Colors = myColors

## would need to order them based on whether they are up or down. Maybe separate them, then patch
## add label for timepoint at which this happens

low_plot <- low_tax_data %>%
    ggplot(aes(y = factor(Genotype, levels = c("BI1", "WT")), x = Abundance, color = Compartment)) +
    geom_boxplot(outlier.shape = NA, width = .5) + # outlier.size = 3, outlier.colour = "black", ) ++
    coord_cartesian(xlim = quantile(low_tax_data$Abundance, c(0, 0.99)))+
    geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
    facet_grid(Taxon ~ Compartment, scales = "free_x") +
    scale_color_manual(values = Colors) +
    theme_bw() +
    theme(
        strip.text = element_text(size = 18), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)
    ) +
    guides(x = guide_axis(angle = 45)) +
    labs(y = "", x = "") 

mid_plot <- mid_tax_data %>%
    ggplot(aes(y = factor(Genotype, levels = c("BI1", "WT")), x = Abundance, color = Compartment)) +
    geom_boxplot(outlier.shape = NA, width = .5) + # outlier.size = 3, outlier.colour = "black", ) ++
    coord_cartesian(xlim = quantile(mid_tax_data$Abundance, c(0, 0.99)))+
    geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
    facet_grid(Taxon ~ Compartment, scales = "free_x") +
    scale_color_manual(values = Colors) +
    theme_bw() +
    theme(
        strip.text = element_text(size = 18), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)
    ) +
    guides(x = guide_axis(angle = 45)) +
    labs(y = "", x = "") 

fusarium_plot <- fusarium_data %>%
    ggplot(aes(y = factor(Genotype, levels = c("BI1", "WT")), x = Abundance, color = Compartment)) +
    geom_boxplot(outlier.shape = NA, width = .5) + # outlier.size = 3, outlier.colour = "black", ) ++
    coord_cartesian(xlim = quantile(fusarium_data$Abundance, c(0, 0.99)))+
    geom_point(aes(shape = Genotype), alpha = .7, size = 3) +
    facet_grid(Taxon ~ Compartment, scales = "free_x") +
    scale_color_manual(values = Colors) +
    theme_bw() +
    theme(
        strip.text.x = element_blank(), strip.text.y = element_text(size = 18), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)
    ) +
    guides(x = guide_axis(angle = 45)) +
    labs(y = "", x = "Relative abundance") 
#ggsave(relabundance, file = "~/Documents/projects/ASP/fungi/data_analysis/root_relative_abundance_significant_differences.pdf", height = 20, width = 15, limitsize = FALSE)

alpha_genus <- alpha_plots[["genus"]]

library(patchwork)
composite_plot <- (low_plot / mid_plot / fusarium_plot + plot_layout(guides = 'collect', heights = c(4,5,1))) | alpha_genus 
composite_plot
ggsave(composite_plot, file = "~/Documents/projects/ASP/fungi/data_analysis/root_relative_abundance_significant_differences.pdf", height = 20, width = 15, limitsize = FALSE)

## maybe i can add one more plot below alpha_genus...




############# Phylogenetic Alpha-diversity #############
library(abdiv)
library(picante) # contains mpd, which I can't find in abdiv
## available alpha-div indices: abdiv::faith_pd(), abdiv::shannon(), abdiv::simpson(), picante::mpd()
## available beta-div indices : abdiv::bray_curtis(), abdiv::unifrac() (unweighted AND weighted), abdiv::phylosor() 
# phylosor is basically the same thing as unifrac, but there's a factor 2 difference
# BUT, for alpha for example, I don't need to compute alpha-phylo for each taxon level, because it uses the phylogenetic tree which is created based on the ASVs.

GetAlphaList <- function(Metadata, FilteredList, TaxonLevels, index) {
    AlphaList <- list()
    for (i in 1:length(TaxonLevels)) {
        level <- TaxonLevels[i]
        data <- rownames_to_column(data.frame(FilteredList[[level]]), var = "Sample_Name")
        index_data <- setNames(data.frame(cbind(data[, 1], diversity(data[, -1], index))), nm = c("Sample_Name", index))
        index_data <- merge(Metadata, index_data, by = "Sample_Name")
        index_data[, ncol(index_data)] <- as.numeric(index_data[, ncol(index_data)])

        kw_test <- index_data %>%
            group_by(Compartment, Timepoint) %>%
            group_map(~ KW_wrapper(.y, .x$shannon, .x$Genotype)) %>%
            bind_rows()

        index_data <- index_data %>%
            dplyr::left_join(kw_test, by = c("Compartment", "Timepoint"))
        AlphaList[[level]] <- index_data
    }
    
    return(AlphaList)
}