# Groundhog
### Abstract


# First run of code = download metadata/files/ load libraries 

```
# Set working directory (optional, but recommended)
setwd("/storage/work/akp5549/groundhog")

#install packages needed
library(phyloseq)
library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
#library(microViz)
library(stringr)
library(tidylog)
library(tidyverse)
library(data.table)
library(writexl)

#### setting paths and such
# ---- Paths ----
RDATA_path <- "/storage/work/akp5549/groundhog/ANSC456-2/ANSC456-2/3.0-dada2/ps_object/by_student"
META_path  <- "/storage/work/akp5549/groundhog/metadata_ANSC456-2.csv"


# Load phyloseq object
# -----------------------------
ps <- readRDS(file.path(RDATA_path, "ps_Groundhog.rds"))


# Load and attach metadata
# -----------------------------
metadata_df <- read.csv(
  META_path,
  header = TRUE,
  stringsAsFactors = FALSE
)

# Check sample names
sample_names(ps)[1:5]
head(metadata_df)

# Make sure rownames of metadata match sample names
rownames(metadata_df) <- metadata_df$SampleID  # change "SampleID" if needed

# Attach metadata to phyloseq object
sample_data(ps) <- sample_data(metadata_df)

```
# error here = 
```
> sample_data(ps) <- sample_data(metadata_df)
Error in validObject(.Object) : invalid class “phyloseq” object: 
 Component sample names do not match.
 Try sample_names()
```


#turns out megtadata was already attached to phylo seq object, great! let's go to next step

#tried to run alpha diversity but error with formula...

```
## -----run diversity metrics on sample content type (area of GI tract)---
unique(sample_data(ps)$Sample_Content)


#------- Alpha diversity -------
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Observed"))

# Add to sample data
sample_data(ps)$Shannon  <- alpha_div$Shannon
sample_data(ps)$Simpson  <- alpha_div$Simpson
sample_data(ps)$Observed <- alpha_div$Observed

library(ggplot2)

ggplot(sample_data(ps), aes(x = Sample_Content, y = Shannon, fill = Sample_Content)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  theme_classic() +
  scale_fill_viridis_d(begin = 0.3, end = 0.85) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  labs(title = "Shannon Diversity Across Groundhog GI Regions",
       x = "Sample Type",
       y = "Shannon Diversity")


### run stats
kruskal.test(Shannon ~ Sample_Content, data = sample_data(ps))

pairwise.wilcox.test(sample_data(ps)$Shannon,
                     sample_data(ps)$Sample_Content,
                     p.adjust.method = "BH")
```

# OOPS we forgot to remove the controls! go back and get that shi out of here
## add before alpha diversity steps:
```
ps <- subset_samples(ps, Sample_Content != "Control 1")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
```


# ok so we still have an error with running stats on alpha diversity...but we did catch that control so good

# we need to convert data frame first!

## so we change to:
```
###-------- run alpha stats ---------------
meta_df <- data.frame(sample_data(ps))

kruskal.test(Shannon ~ Sample_Content, data = meta_df)

pairwise.wilcox.test(sample_data(ps)$Shannon,
                     sample_data(ps)$Sample_Content,
                     p.adjust.method = "BH")
```


# trouble with microviz so we changed beta diversity section to:
```
# ------------ Beta diversity ------------
# Bray-Curtis distance
bray_dist <- phyloseq::distance(ps, method = "bray")

# PCoA
pcoa_res <- ordinate(ps, method = "PCoA", distance = bray_dist)

# Plot
plot_ordination(ps, pcoa_res, color = "Sample_Content") +
  geom_point(size = 4, alpha = 0.7) +
  theme_classic() +
  scale_color_viridis_d() +
  labs(title = "Bray-Curtis PCoA by Groundhog GI Region")
```


# but then we dont have ellipses over points, so we change to
```
# Bray-Curtis distance and PCoA
bray_dist <- phyloseq::distance(ps, method = "bray")
pcoa_res <- ordinate(ps, method = "PCoA", distance = bray_dist)

# Convert sample_data to dataframe
meta_df <- as.data.frame(sample_data(ps))

# Plot with ellipses
plot_ordination(ps, pcoa_res, type = "samples", color = "Sample_Content") +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(aes(group = Sample_Content, color = Sample_Content), linetype = 2, level = 0.95) +
  theme_classic() +
  scale_color_viridis_d() +
  labs(title = "Bray-Curtis PCoA with Ellipses by GI Region")
```

# it works! Now onto relative abundance 

# Here is phyla and genus level relative abundance plots and such
```


# -------------- relative abundance -------------------------

# Function to calculate Top-N taxa per Sample_Content
# =============================
topN_relative_abundance <- function(ps, taxrank = "Genus", topN = 10) {
  
  # Collapse to taxonomic rank
  ps.rank <- tax_glom(ps, taxrank = taxrank, NArm = FALSE)
  
  # Remove samples with zero counts
  ps.rank <- prune_samples(sample_sums(ps.rank) > 0, ps.rank)
  
  # Transform to relative abundance
  ps.rel <- transform_sample_counts(ps.rank, function(x) x / sum(x))
  df <- psmelt(ps.rel)
  
  # Clean taxonomy
  df[[taxrank]] <- as.character(df[[taxrank]])
  df[[taxrank]][is.na(df[[taxrank]]) |
                  grepl("^\\[.*\\]|Candidatus|Incertae Sedis|UCG|group|^[A-Z0-9_-]+$", df[[taxrank]])] <- "Unclassified"
  
  # Loop over Sample_Contents separately
  out <- lapply(unique(df$Sample_Content), function(g) {
    
    df_g <- df  filter(Sample_Content == g)
    
    # Top N taxa per group (excluding Unclassified)
    top_taxa <- df_g 
      filter(.data[[taxrank]] != "Unclassified") 
      group_by(.data[[taxrank]]) 
      summarise(Total = sum(Abundance), .groups = "drop") 
      slice_max(Total, n = topN) 
      pull(.data[[taxrank]])
    
    # Assign Other / Unclassified
    df_g <- df_g 
      mutate(
        Taxon_plot = case_when(
          .data[[taxrank]] %in% top_taxa ~ .data[[taxrank]],
          .data[[taxrank]] == "Unclassified" ~ "Unclassified",
          TRUE ~ "Other"
        )
      )
    
    # Sum abundances per sample & enforce factor order
    df_g <- df_g 
      group_by(Sample, Sample_Content, Taxon_plot) 
      summarise(Abundance = sum(Abundance), .groups = "drop") 
      mutate(Taxon_plot = factor(Taxon_plot, levels = c(sort(top_taxa), "Other", "Unclassified")))
    
    df_g
  })
  
  # Combine all groups
  bind_rows(out)
}

# =============================
# Run Top 10 analyses
# =============================
genus_rel  <- topN_relative_abundance(ps, taxrank = "Genus", topN = 10)
phylum_rel <- topN_relative_abundance(ps, taxrank = "Phylum", topN = 10)


# -----------------------------
# Create consistent color palette
# =============================
create_taxon_colors <- function(df) {
  all_taxa <- unique(df$Taxon_plot)
  top_taxa <- setdiff(all_taxa, c("Other", "Unclassified"))
  
  # Use qualitative palette from RColorBrewer
  n_top <- length(top_taxa)
  if(n_top <= 12){
    top_colors <- brewer.pal(n_top, "Set3")
  } else {
    top_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_top)
  }
  
  taxon_colors <- c(
    setNames(top_colors, top_taxa),
    "Other" = "gray50",
    "Unclassified" = "yellow"
  )
  
  # Ensure factor levels match colors
  df$Taxon_plot <- factor(df$Taxon_plot, levels = names(taxon_colors))
  
  list(df = df, taxon_colors = taxon_colors)
}

genus_pal <- create_taxon_colors(genus_rel)
genus_rel <- genus_pal$df
genus_colors <- genus_pal$taxon_colors

phylum_pal <- create_taxon_colors(phylum_rel)
phylum_rel <- phylum_pal$df
phylum_colors <- phylum_pal$taxon_colors

# =============================
# Function to plot per Sample_Content
# =============================
plot_topN <- function(df, taxname = "Genus", colors) {
  for (g in unique(df$Sample_Content)) {
    
    df_g <- df  filter(Sample_Content == g)
    
    p <- ggplot(df_g, aes(x = Sample, y = Abundance, fill = Taxon_plot)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = colors) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid.major.x = element_blank()) +
      labs(title = paste("Top 10", taxname, "-", g),
           y = "Relative Abundance",
           fill = taxname)
    
    print(p)
  }
}

# =============================
# Generate plots
# =============================
plot_topN(genus_rel, "Genus", genus_colors)
plot_topN(phylum_rel, "Phylum", phylum_colors)

# =============================
# Save PDF of all plots
# =============================
pdf("Top10_relative_abundance_plots.pdf", width = 12, height = 8)
plot_topN(genus_rel, "Genus", genus_colors)
plot_topN(phylum_rel, "Phylum", phylum_colors)
dev.off()


# =============================
# Save summary table (mean % + SD) by Sample_Content 
# =============================
genus_summary <- genus_rel 
  group_by(Sample_Content, Taxon_plot) 
  summarise(MeanPercent = mean(Abundance)*100,
            SD = sd(Abundance)*100,
            .groups = "drop") 
  arrange(Sample_Content, desc(MeanPercent)) 
  mutate(MeanPercent = round(MeanPercent,2),
         SD = round(SD,2))

write.csv(genus_summary, "Top10_Genus_mean_SD.csv", row.names = FALSE)

# Optional: save as Excel
writexl::write_xlsx(list(Genus_summary = genus_summary), "Top10_Genus_mean_SD.xlsx")
```
