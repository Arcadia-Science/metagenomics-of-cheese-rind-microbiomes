library(tidyverse)
library(sourmashconsumr)
library(ggpubr)

#########################################################
# TMI Illumina sourmashconsumr exploration
#########################################################

## read in reads files for TMI illumina samples 
# reads signatures
tmi_illumina_sigs_directory <- ("processed_data/2023-05-09-TMI-processed-data/TMI_illumina/sketch/")
tmi_illumina_sigs_files <- list.files(path = tmi_illumina_sigs_directory, pattern = "*.reads.sig",  full.names = TRUE)
tmi_illumina_reads_sigs <- read_signature(tmi_illumina_sigs_files)

# reads compare
tmi_illumina_reads_compare_csv <- read_compare_csv("processed_data/2023-05-09-TMI-processed-data/TMI_illumina/compare/reads.comp.csv", sample_to_rownames = F)

# reads gather CSVs
tmi_illumina_gather_directory <- ("processed_data/2023-05-09-TMI-processed-data/TMI_illumina/gather/")
tmi_illumina_gather_csvs <- list.files(path = tmi_illumina_gather_directory, pattern = "*.reads.gather.csv", full.names = TRUE)
tmi_illumina_reads_gather <- read_gather(tmi_illumina_gather_csvs, intersect_bp_threshold = 50000)

# reads taxonomy annotate
tmi_illumina_taxonomy_csvs <- list.files(path = "processed_data/2023-05-09-TMI-processed-data/TMI_illumina/taxonomy/", pattern = "*.reads.gather.with-lineages.csv.gz", full.names = TRUE)
tmi_illumina_reads_taxonomy <- read_taxonomy_annotate(tmi_illumina_taxonomy_csvs, intersect_bp_threshold = 50000, separate_lineage = T)

## read in assemblies files for TMI illumina samples 
# assembly signatures
tmi_illumina_assembs_sigs_files <- list.files(path = tmi_illumina_sigs_directory, pattern = "*.assembly.sig",  full.names = TRUE)
tmi_illumina_assembs_sigs <- read_signature(tmi_illumina_assembs_sigs_files)

# assemblies compare
tmi_illumina_assemblies_compare_csv <- read_compare_csv("processed_data/2023-05-09-TMI-processed-data/TMI_illumina/compare/assembly.comp.csv", sample_to_rownames = F)

# assemblies gather
tmi_illumina_assembs_gather_csvs <- list.files(path = tmi_illumina_gather_directory, pattern = "*.assembly.gather.csv", full.names = TRUE)
tmi_illumina_assembs_gather <- read_gather(tmi_illumina_assembs_gather_csvs, intersect_bp_threshold = 50000)

# assemblies taxonomy 
tmi_illumina_assembs_taxonomy_csvs <- list.files(path = "processed_data/2023-05-09-TMI-processed-data/TMI_illumina/taxonomy/", pattern = "*.assembly.gather.with-lineages.csv.gz", full.names = TRUE)
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
illumina_assemb_classified_plot <- tmi_illumina_assembs_gather %>% 
  filter(!grepl("plant", filename)) %>% 
  plot_gather_classified() +
  ggtitle("Classified Sequences in Illumina Assemblies")


# plotting taxonomy annotate results
plot_taxonomy_annotate_sankey(tmi_illumina_reads_taxonomy, tax_glom_level = "order")
plot_taxonomy_annotate_sankey(tmi_illumina_assembs_taxonomy, tax_glom_level = "order")

#########################################################
# TMI Nanopore sourmashconsumr exploration
#########################################################

## read in reads files for TMI Nanopore samples 
# reads signatures
tmi_nanopore_sigs_directory <- ("processed_data/2023-05-09-TMI-processed-data/TMI_nanopore/sketch/")
tmi_nanopore_sigs_files <- list.files(path = tmi_nanopore_sigs_directory, pattern = "*.reads.sig",  full.names = TRUE)
tmi_nanopore_reads_sigs <- read_signature(tmi_nanopore_sigs_files)

# reads compare
tmi_nanopore_reads_compare_csv <- read_compare_csv("processed_data/2023-05-09-TMI-processed-data/TMI_nanopore/compare/reads.comp.csv", sample_to_rownames = F)

# reads gather CSVs
tmi_nanopore_gather_directory <- ("processed_data/2023-05-09-TMI-processed-data/TMI_nanopore/gather/")
tmi_nanopore_gather_csvs <- list.files(path = tmi_nanopore_gather_directory, pattern = "*.reads.gather.csv", full.names = TRUE)
tmi_nanopore_reads_gather <- read_gather(tmi_nanopore_gather_csvs, intersect_bp_threshold = 50000)

# reads taxonomy annotate
tmi_nanopore_taxonomy_csvs <- list.files(path = "processed_data/2023-05-09-TMI-processed-data/TMI_nanopore/taxonomy/", pattern = "*.reads.gather.with-lineages.csv.gz", full.names = TRUE)
tmi_nanopore_reads_taxonomy <- read_taxonomy_annotate(tmi_nanopore_taxonomy_csvs, intersect_bp_threshold = 50000, separate_lineage = T)

## read in assemblies files for TMI Nanopore samples 
# assembly signatures
tmi_nanopore_assembs_sigs_files <- list.files(path = tmi_nanopore_sigs_directory, pattern = "*.assembly.sig",  full.names = TRUE)
tmi_nanopore_assembs_sigs <- read_signature(tmi_nanopore_assembs_sigs_files)

# assemblies compare
tmi_nanopore_assemblies_compare_csv <- read_compare_csv("processed_data/2023-05-09-TMI-processed-data/TMI_nanopore/compare/assembly.comp.csv", sample_to_rownames = F)

# assemblies gather
tmi_nanopore_assembs_gather_csvs <- list.files(path = tmi_nanopore_gather_directory, pattern = "*.assembly.gather.csv", full.names = TRUE)
tmi_nanopore_assembs_gather <- read_gather(tmi_nanopore_assembs_gather_csvs, intersect_bp_threshold = 50000)

# assemblies taxonomy 
tmi_nanopore_assembs_taxonomy_csvs <- list.files(path = "processed_data/2023-05-09-TMI-processed-data/TMI_nanopore/taxonomy/", pattern = "*.assembly.gather.with-lineages.csv.gz", full.names = TRUE)
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
nanopore_assembs_classified_plot <- plot_gather_classified(tmi_nanopore_assembs_gather) +
  ggtitle("Classified Sequences in Nanopore Assemblies")

# plotting taxonomy annotate results
plot_taxonomy_annotate_sankey(tmi_nanopore_reads_taxonomy, tax_glom_level = "order")
nanopore_assembs_tax_sankey_plot <- plot_taxonomy_annotate_sankey(tmi_nanopore_assembs_taxonomy, tax_glom_level = "order")

# eligo nanopore time series plots 
eligo_taxonomy <- tmi_nanopore_assembs_taxonomy %>% 
  filter(grepl('Eligo', query_name))
eligo_time_df <- data.frame(query_name = c('Eligo2weeks_T1.assembly', 'Eligo4weeks_T1.assembly', 'Eligo12weeks_T1.assembly'),
                            time = c('1', '2', '3'))
plot_taxonomy_annotate_ts_alluvial(eligo_taxonomy, eligo_time_df, tax_glom_level = "phylum", fraction_threshold = 0.01) +
  ggplot2::scale_fill_brewer(palette = "Paired")

# whitney nanopore time series plots
whitney_taxonomy <- tmi_nanopore_assembs_taxonomy %>% 
  filter(grepl('Whitney', query_name))
whitney_time_df <- data.frame(query_name = c('Whitney1month_T1.assembly', 'Whitney2months_T1.assembly', 'Whitney4months_T1.assembly'),
                              time = c('1', '2', '3'))
plot_taxonomy_annotate_ts_alluvial(whitney_taxonomy, whitney_time_df, tax_glom_level="phylum", fraction_threshold = 0.01)

# oma nanopore time series plots
oma_taxonomy <- tmi_nanopore_assembs_taxonomy %>% 
  filter(grepl('Oma', query_name))
oma_time_df <- data.frame(query_name = c('Oma2weeks_T1.assembly', 'Oma4weeks_T1.assembly', 'Oma8weeks_T1.assembly'),
                          time = c('1', '2', '3'))
plot_taxonomy_annotate_ts_alluvial(oma_taxonomy, oma_time_df, tax_glom_level = "phylum", fraction_threshold = 0.01)

#########################################################
# Plot grids and save figures
#########################################################


classified_plot <- ggarrange(illumina_assemb_classified_plot, nanopore_assembs_classified_plot, ncol=2, labels = c("A", "B"), common.legend = TRUE, legend = "bottom", align = "h")

ggsave("figures/TMI_classified_plot.png", classified_plot, width=30, height=15, units=c("cm"))

ggsave("figures/TMI_nanopore_assembs_tax_sankey.png", nanopore_assembs_tax_sankey_plot, width=30, height=20, units=c("cm"))
