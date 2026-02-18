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
