library(tidyverse)
library(sourmashconsumr)
library(ggpubr)
library(ArcadiaColorBrewer)

#########################################################
# TMI Illumina sourmashconsumr exploration
#########################################################

## read in reads files for TMI illumina samples 
# reads signatures
tmi_illumina_sigs_directory <- ("processed_data/2023-05-23-TMI-processed-data/TMI_illumina/sketch/")
tmi_illumina_sigs_files <- list.files(path = tmi_illumina_sigs_directory, pattern = "*.reads.sig",  full.names = TRUE)
tmi_illumina_reads_sigs <- read_signature(tmi_illumina_sigs_files)

# reads compare
tmi_illumina_reads_compare_csv <- read_compare_csv("processed_data/2023-05-23-TMI-processed-data/TMI_illumina/compare/reads.comp.csv", sample_to_rownames = F)

# reads gather CSVs
tmi_illumina_gather_directory <- ("processed_data/2023-05-23-TMI-processed-data/TMI_illumina/gather/")
tmi_illumina_gather_csvs <- list.files(path = tmi_illumina_gather_directory, pattern = "*.reads.gather.csv", full.names = TRUE)
tmi_illumina_reads_gather <- read_gather(tmi_illumina_gather_csvs, intersect_bp_threshold = 50000)

# reads taxonomy annotate
tmi_illumina_taxonomy_csvs <- list.files(path = "processed_data/2023-05-23-TMI-processed-data/TMI_illumina/taxonomy/", pattern = "*.reads.gather.with-lineages.csv.gz", full.names = TRUE)
tmi_illumina_reads_taxonomy <- read_taxonomy_annotate(tmi_illumina_taxonomy_csvs, intersect_bp_threshold = 50000, separate_lineage = T)

## read in assemblies files for TMI illumina samples 
# assembly signatures
tmi_illumina_assembs_sigs_files <- list.files(path = tmi_illumina_sigs_directory, pattern = "*.assembly.sig",  full.names = TRUE)
tmi_illumina_assembs_sigs <- read_signature(tmi_illumina_assembs_sigs_files)

# assemblies compare
tmi_illumina_assemblies_compare_csv <- read_compare_csv("processed_data/2023-05-23-TMI-processed-data/TMI_illumina/compare/assembly.comp.csv", sample_to_rownames = F)

# assemblies gather
tmi_illumina_assembs_gather_csvs <- list.files(path = tmi_illumina_gather_directory, pattern = "*.assembly.gather.csv", full.names = TRUE)
tmi_illumina_assembs_gather <- read_gather(tmi_illumina_assembs_gather_csvs, intersect_bp_threshold = 50000)

# assemblies taxonomy 
tmi_illumina_assembs_taxonomy_csvs <- list.files(path = "processed_data/2023-05-23-TMI-processed-data/TMI_illumina/taxonomy/", pattern = "*.assembly.gather.with-lineages.csv.gz", full.names = TRUE)
tmi_illumina_assembs_taxonomy <- read_taxonomy_annotate(tmi_illumina_assembs_taxonomy_csvs, intersect_bp_threshold = 50000, separate_lineage = T)

## plotting TMI illumina reads and assemblies results
# plot compare mds and heatmap
tmi_illumina_compare_reads_mds_df <- make_compare_mds(tmi_illumina_reads_compare_csv)
tmi_illumina_compare_assemblies_mds_df <- make_compare_mds(tmi_illumina_assemblies_compare_csv)

plot_compare_mds(tmi_illumina_compare_reads_mds_df)
plot_compare_mds(tmi_illumina_compare_assemblies_mds_df)
plot_compare_heatmap(tmi_illumina_reads_compare_csv, cexRow = 0.75, cexCol = 0.75)
illumina_assemblies_compared_plot <- plot_compare_heatmap(tmi_illumina_assemblies_compare_csv, cexRow = 0.75, cexCol = 0.75)

# plotting gather results
plot_gather_classified(tmi_illumina_reads_gather)

plot_gather_classified(tmi_illumina_assembs_gather)

# remove low hits to plant db
arcadia.pal(6, 'Accent')
illumina_assemb_classified_plot <- tmi_illumina_assembs_gather %>% mutate(query_name = gsub(".assembly", "", query_name)) %>% 
  plot_gather_classified() +
  scale_fill_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846")) +
  ggtitle("Classified Sequences in Illumina Assemblies") + 
  scale_y_continuous(expand = c(0,0)) 

illumina_assemb_classified_plot


# plotting taxonomy annotate results
plot_taxonomy_annotate_sankey(tmi_illumina_reads_taxonomy, tax_glom_level = "order")
plot_taxonomy_annotate_sankey(tmi_illumina_assembs_taxonomy, tax_glom_level = "order")

#########################################################
# TMI Nanopore sourmashconsumr exploration
#########################################################

## read in reads files for TMI Nanopore samples 
# reads signatures
tmi_nanopore_sigs_directory <- ("processed_data/2023-05-23-TMI-processed-data/TMI_nanopore/sketch/")
tmi_nanopore_sigs_files <- list.files(path = tmi_nanopore_sigs_directory, pattern = "*.reads.sig",  full.names = TRUE)
tmi_nanopore_reads_sigs <- read_signature(tmi_nanopore_sigs_files)

# reads compare
tmi_nanopore_reads_compare_csv <- read_compare_csv("processed_data/2023-05-23-TMI-processed-data/TMI_nanopore/compare/reads.comp.csv", sample_to_rownames = F)

# reads gather CSVs
tmi_nanopore_gather_directory <- ("processed_data/2023-05-23-TMI-processed-data/TMI_nanopore/gather/")
tmi_nanopore_gather_csvs <- list.files(path = tmi_nanopore_gather_directory, pattern = "*.reads.gather.csv", full.names = TRUE)
tmi_nanopore_reads_gather <- read_gather(tmi_nanopore_gather_csvs, intersect_bp_threshold = 50000)

# reads taxonomy annotate
tmi_nanopore_taxonomy_csvs <- list.files(path = "processed_data/2023-05-23-TMI-processed-data/TMI_nanopore/taxonomy/", pattern = "*.reads.gather.with-lineages.csv.gz", full.names = TRUE)
tmi_nanopore_reads_taxonomy <- read_taxonomy_annotate(tmi_nanopore_taxonomy_csvs, intersect_bp_threshold = 50000, separate_lineage = T)

## read in assemblies files for TMI Nanopore samples 
# assembly signatures
tmi_nanopore_assembs_sigs_files <- list.files(path = tmi_nanopore_sigs_directory, pattern = "*.assembly.sig",  full.names = TRUE)
tmi_nanopore_assembs_sigs <- read_signature(tmi_nanopore_assembs_sigs_files)

# assemblies compare
tmi_nanopore_assemblies_compare_csv <- read_compare_csv("processed_data/2023-05-23-TMI-processed-data/TMI_nanopore/compare/assembly.comp.csv", sample_to_rownames = F)

# assemblies gather
tmi_nanopore_assembs_gather_csvs <- list.files(path = tmi_nanopore_gather_directory, pattern = "*.assembly.gather.csv", full.names = TRUE)
tmi_nanopore_assembs_gather <- read_gather(tmi_nanopore_assembs_gather_csvs, intersect_bp_threshold = 50000)

# assemblies taxonomy 
tmi_nanopore_assembs_taxonomy_csvs <- list.files(path = "processed_data/2023-05-23-TMI-processed-data/TMI_nanopore/taxonomy/", pattern = "*.assembly.gather.with-lineages.csv.gz", full.names = TRUE)
tmi_nanopore_assembs_taxonomy <- read_taxonomy_annotate(tmi_nanopore_assembs_taxonomy_csvs, intersect_bp_threshold = 50000, separate_lineage = T)

## plotting TMI Nanopore reads and assemblies results
# plot compare mds and heatmap
tmi_nanopore_compare_reads_mds_df <- make_compare_mds(tmi_nanopore_reads_compare_csv)
tmi_nanopore_compare_assemblies_mds_df <- make_compare_mds(tmi_nanopore_assemblies_compare_csv)

plot_compare_mds(tmi_nanopore_compare_reads_mds_df)
plot_compare_mds(tmi_nanopore_compare_assemblies_mds_df)
plot_compare_heatmap(tmi_nanopore_reads_compare_csv, cexRow = 0.75, cexCol = 0.75)
plot_compare_heatmap(tmi_nanopore_assemblies_compare_csv, cexRow = 0.75, cexCol = 0.75)

# plotting gather results
plot_gather_classified(tmi_nanopore_reads_gather)
nanopore_assembs_classified_plot <- tmi_nanopore_assembs_gather %>%  mutate(query_name = gsub(".assembly", "", query_name)) %>% 
  plot_gather_classified() +
  scale_fill_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846")) +
  ggtitle("Classified Sequences in Nanopore Assemblies") + 
  scale_y_continuous(expand = c(0,0)) 

nanopore_assembs_classified_plot

# plotting taxonomy annotate results
options(repr.plot.width = 8.5, repr.plot.height = 3, repr.plot.res = 300)

nanopore_assembs_tax_sankey_plot <- plot_taxonomy_annotate_sankey(tmi_nanopore_assembs_taxonomy, tax_glom_level = "order", palette = grDevices::colorRampPalette(c("#C6E7F4", "#F8C5C1", "#F5E4BE", "#B5BEA4", 
                                                                                                                                         "#DCBFFC", "#B6C8D4", "#DA9085",
                                                                                                                                         "#F5CBE4", "#BABEE0", "#D1EADF", "#F1E8DA"))(n = 11),
                                                                  label = F) +
  ggforce::geom_parallel_sets_labels(colour = 'black', angle = 360, size = 3, fontface = "italic", hjust = -0.25) +
  labs(x = "Taxonomic rank") +
  scale_x_continuous(labels = c("Domain", "Phylum", "Class", "Order", ""),
                     breaks = c(1, 2, 3, 4, 5),
                     limits = c(.75, 5)) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 16))

nanopore_assembs_tax_sankey_plot

# eligo nanopore time series plots 
eligo_taxonomy <- tmi_nanopore_assembs_taxonomy %>% 
  filter(grepl('EL', query_name))
eligo_time_df <- data.frame(query_name = c('EL2weeks.assembly', 'EL4weeks.assembly', 'EL12weeks.assembly'),
                            time = c('1', '2', '3'))
plot_taxonomy_annotate_ts_alluvial(eligo_taxonomy, eligo_time_df, tax_glom_level = "order", fraction_threshold = 0.01) +
  ggplot2::scale_fill_brewer(palette = "Paired")

#########################################################
# Plot grids and save figures
#########################################################


classified_plot <- ggarrange(illumina_assemb_classified_plot, nanopore_assembs_classified_plot, ncol=2, labels = c("A", "B"), common.legend = TRUE, legend = "bottom", align = "h")

classified_plot

ggsave("figures/TMI_classified_plot.png", classified_plot, width=30, height=15, units=c("cm"))

ggsave("figures/TMI_nanopore_assembs_tax_sankey.png", nanopore_assembs_tax_sankey_plot, width=30, height=20, units=c("cm"))
