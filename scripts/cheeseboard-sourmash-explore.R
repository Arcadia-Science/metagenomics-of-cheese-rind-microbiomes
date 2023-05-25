library(tidyverse)
library(sourmashconsumr)

#########################################################
# Sourmash results on cheeseboard reads 
#########################################################

## read in all files for cheeseboard reads

# reads signatures
cheeseboard_sigs_directory <- ("processed_data/2023-04-24-cheeseboard-processed-data/sourmash_files/signatures")
cheeseboard_sigs_files <- list.files(path = cheeseboard_sigs_directory, pattern = "*.reads.sig",  full.names = TRUE)
cheeseboard_reads_sigs <- read_signature(cheeseboard_sigs_files)

# reads compare
cheeseboard_reads_compare_csv <- read_compare_csv("processed_data/2023-04-24-cheeseboard-processed-data/sourmash_files/compare/reads.comp.csv", sample_to_rownames = F)

# reads gather CSVs
cheeseboard_gather_directory <- ("processed_data/2023-04-24-cheeseboard-processed-data/sourmash_files/gather")
cheeseboard_gather_csvs <- list.files(path = cheeseboard_gather_directory, pattern = "*.reads.gather.csv", full.names = TRUE)
cheeseboard_reads_gather <- read_gather(cheeseboard_gather_csvs, intersect_bp_threshold = 50000)

# reads taxonomy annotate
cheeseboard_taxonomy_csvs <- list.files(path = cheeseboard_gather_directory, pattern = "*.reads.gather.with-lineages.csv.gz", full.names = TRUE)
cheeseboard_reads_taxonomy <- read_taxonomy_annotate(cheeseboard_taxonomy_csvs, intersect_bp_threshold = 50000, separate_lineage = T)

# plot compare mds and heatmap
cheeseboard_compare_reads_mds_df <- make_compare_mds(cheeseboard_reads_compare_csv)
plot_compare_mds(cheeseboard_compare_reads_mds_df)
plot_compare_heatmap(cheeseboard_reads_compare_csv, cexRow = 0.75, cexCol = 0.75)

# plotting gather results
plot_gather_classified(cheeseboard_reads_gather)

# plotting taxonomy annotate results
glom_taxonomy_reads <- tax_glom_taxonomy_annotate(cheeseboard_reads_taxonomy, tax_glom_level = "genus", glom_var = "f_unique_to_query")
glom_taxonomy_reads_upset <- from_taxonomy_annotate_to_upset_inputs(cheeseboard_reads_taxonomy, tax_glom_level = "order")
plot_taxonomy_annotate_upset(glom_taxonomy_reads_upset, fill = "phylum")
plot_taxonomy_annotate_sankey(cheeseboard_reads_taxonomy, tax_glom_level = "order")


#########################################################
# Sourmash results on cheeseboard assemblies
#########################################################

## read in all files for cheeseboard assemblies

# assemblies signatures
cheeseboard_sigs_assembs_files <- list.files(path = cheeseboard_sigs_directory, pattern = "*.assembly.sig",  full.names = TRUE)
cheeseboard_assemblies_sigs <- read_signature(cheeseboard_sigs_assembs_files)

# assemblies compare
cheeseboard_assemblies_compare_csv <- read_compare_csv("processed_data/2023-04-24-cheeseboard-processed-data/sourmash_files/compare/assembly.comp.csv", sample_to_rownames = F)

# assemblies gather CSVs
cheeseboard_gather_assemblies_csvs <- list.files(path = cheeseboard_gather_directory, pattern = "*.assembly.gather.csv", full.names = TRUE)
cheeseboard_assemblies_gather <- read_gather(cheeseboard_gather_assemblies_csvs, intersect_bp_threshold = 50000)

# assemblies taxonomy annotate
cheeseboard_assemblies_taxonomy_csvs <- list.files(path = cheeseboard_gather_directory, pattern = "*.assembly.gather.with-lineages.csv.gz", full.names = TRUE)
cheeseboard_assemblies_taxonomy <- read_taxonomy_annotate(cheeseboard_assemblies_taxonomy_csvs, intersect_bp_threshold = 50000, separate_lineage = T)

# plot compare mds and heatmap
cheeseboard_compare_assemblies_mds_df <- make_compare_mds(cheeseboard_assemblies_compare_csv)
plot_compare_mds(cheeseboard_compare_assemblies_mds_df)
plot_compare_heatmap(cheeseboard_assemblies_compare_csv, cexRow = 0.75, cexCol = 0.75)

# plotting gather results
plot_gather_classified(cheeseboard_assemblies_gather)

# plotting taxonomy annotate results
glom_taxonomy_assemblies <- tax_glom_taxonomy_annotate(cheeseboard_assemblies_taxonomy, tax_glom_level = "genus", glom_var = "f_unique_to_query")
glom_taxonomy_assemblies_upset <- from_taxonomy_annotate_to_upset_inputs(cheeseboard_assemblies_taxonomy, tax_glom_level = "order")
plot_taxonomy_annotate_sankey(cheeseboard_assemblies_taxonomy, tax_glom_level = "order")

