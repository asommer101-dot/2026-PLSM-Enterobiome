###Author: Andrew J Sommer
###Date Modified: 3/17/2026
###Description: Code for Anaylsis of plasmid Abundance

library(dplyr)
library(vegan)
library(tibble) 
library(ggplot2)
library(MetBrewer)
library(cowplot)
library(forcats)
library(cowplot)
library(patchwork)
library(viridis)
library(rstatix)

#####Read in and Merge Data
#Metadata
sample_meta <- read.csv("Import/Metdata/Sample-Metadata.csv")
head(sample_meta)


###############################################Extended Figure 2A (Plasmid Barplot)#############################################

#Read in Data
Type_Abundance <- read.csv("Import/Data/Taxonomic_Abundance/Family_Plasmid_Abundance_50percentcomplete_5percentcontam.csv") 
Type_Abundance <- Type_Abundance %>% left_join(sample_meta, by = "alias") #add metadata
head(Type_Abundance)

unique(Type_Abundance$Taxa)


#### Recode
Type_Abundance <- Type_Abundance %>%
  mutate(
    Taxa = case_when(
      Taxa %in% c("Unassigned", "Unclassified", "Unclassified_Bacteria") ~ "Unclassified",
      Taxa %in% c("Enterobacteriaceae/Unclassified_Bacteria", "Enterobacteriaceae") ~ "Enterobacteriaceae",
      Taxa %in% c("Lachnospiraceae/Unclassified_Bacteria", "Lachnospiraceae") ~ "Lachnospiraceae",
      Taxa %in% c("Bacteroidaceae/Unclassified_Bacteria", "Bacteroidaceae") ~ "Bacteroidaceae",
      Taxa %in% c("Bacteroidaceae/Unclassified_Bacteria", "Bacteroidaceae") ~ "Bacteroidaceae",
      Taxa %in% c("Lactobacillaceae/Veillonellaceae", 
                  "Bacteroidaceae/Desulfovibrionaceae",
                  "Bacteroidaceae/Tannerellaceae",
                  "Acutalibacteraceae/Lachnospiraceae",
                  "Bacteroidaceae/Rikenellaceae/Tannerellaceae",
                  "Desulfovibrionaceae/Lactobacillaceae",
                  "Acutalibacteraceae/Ruminococcaceae",
                  "Bacteroidaceae/Erysipelotrichaceae/Lachnospiraceae",
                  ""
      ) ~ "Mixed",
      TRUE ~ Taxa))
unique(Type_Abundance$Taxa)


Type_Abundance_collapsed <- Type_Abundance %>%
  group_by(alias, Taxa, treatment, ID) %>%
  summarise(
    Read_Count = sum(Read_Count, na.rm = TRUE),
    TPM        = sum(TPM, na.rm = TRUE),
    .groups = "drop"
  )



#####Deterime Relative abundance by TPM and Read Count
Type_Abundance_collapsed <-
  Type_Abundance_collapsed %>%
  group_by(ID) %>%                          # group by sample
  mutate(Count_Relative_Abundance = Read_Count / sum(Read_Count)) %>%  # divide by total per sample
  mutate(TPM_Relative_Abundance = TPM / sum(TPM)) %>%  # divide by total per sample
  ungroup()


######Identify Low Abundance Taxa
Type_mean_abundance <- Type_Abundance_collapsed %>%
  group_by(Taxa) %>%
  summarise(mean_Count_Relative_Abundance = mean(Count_Relative_Abundance, na.rm = TRUE),
            mean_TPM_Relative_Abundance = mean(TPM_Relative_Abundance, na.rm = TRUE)) %>%
  ungroup()

threshold <- 0.15
low_taxa <- Type_mean_abundance  %>%
  filter(mean_Count_Relative_Abundance < threshold) %>%
  pull(Taxa)
low_taxa



#######Filter out low abundance taxa
Type_Abundance_collapsed <- Type_Abundance_collapsed %>%
  mutate(taxa_grouped = ifelse(Taxa %in% low_taxa, "Other", Taxa))

Type_Abundance_lumped <- Type_Abundance_collapsed %>%
  # Lump all “Other” together per sample
  group_by(ID, treatment, taxa_grouped)%>%
  summarise(Count_Relative_Abundance = sum(Count_Relative_Abundance, na.rm = TRUE),
            TPM_Relative_Abundance = sum(TPM_Relative_Abundance, na.rm = TRUE))  %>%
  ungroup()



unique(Type_Abundance_lumped$taxa_grouped)


#####Set Colors and order
Taxa_colors <- c(
  "Enterobacteriaceae" = "#CA483A",
  "Bacteroidaceae"         = "#7B9671",
  "Tannerellaceae"         = "#599675",
  "Rikenellaceae"         = "#969571",
  "Lactobacillaceae"   = "#4A90E2",
  "Enterococcaceae"    = "#ff7f00", 
  "Other"               = "gray60",
  "Mixed"       = "grey80"
)



Type_Abundance_lumped$taxa_grouped <- factor(Type_Abundance_lumped$taxa_grouped, levels = names(Taxa_colors))
unique(Type_Abundance_lumped$taxa_grouped)


Type_Abundance_lumped$ID <- factor(
  Type_Abundance_lumped$ID,
  levels = c(
    sort(grep("^HC", unique(Type_Abundance_lumped$ID), value = TRUE)),
    sort(grep("^LC", unique(Type_Abundance_lumped$ID), value = TRUE)),
    sort(grep("^CD", unique(Type_Abundance_lumped$ID), value = TRUE))
  )
)

Extended_Fig2A <- ggplot(Type_Abundance_lumped,
                               aes(x = ID, y = Count_Relative_Abundance, fill = taxa_grouped)) +
  geom_bar(stat = "identity", color = "black", width = 1, linewidth = 0.05) +  
  ylab("Relative Abundance") +
  xlab(" ") +
  theme_bw() +
  scale_fill_manual(values = Taxa_colors) +
  theme_half_open() +
  scale_x_discrete(expand = c(0, 0)) +  # remove axis padding
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(fill = "Class") +
  theme(
    axis.text.x = element_text(hjust = 1, size = 8, angle = 60, vjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.line = element_line(size = 0.3)
  ) +
  theme(legend.position = "bottom")
Extended_Fig2A


#####Summarize by Avg across samples

Avg_Abundance <- Type_Abundance_collapsed %>%
  group_by(treatment, Taxa) %>%
  summarise(
    TPM_Avg_RelAb = mean(Count_Relative_Abundance, na.rm = TRUE),
    Count_Avg_RelAb = mean(TPM_Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  )


Avg_Abundance <- Avg_Abundance %>%
  mutate(TPM_Avg_RelAb_percent = TPM_Avg_RelAb * 100) %>%
  mutate(Count_Avg_RelAb_percent = Count_Avg_RelAb * 100)

Avg_Abundance


Avg_Abundance_combined <- Type_Abundance_collapsed %>%
  group_by(Taxa) %>%
  summarise(
    TPM_Avg_RelAb = mean(Count_Relative_Abundance, na.rm = TRUE),
    Count_Avg_RelAb = mean(TPM_Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  )
Avg_Abundance_combined

#### Calculate Entero Prev
entero_counts <- Type_Abundance_collapsed %>%
  filter(Taxa == "Enterobacterales" & Read_Count > 0) %>%
  group_by(treatment) %>%
  summarise(n_present = n_distinct(ID))

# Get total samples per treatment from sample_meta
sample_sizes <- sample_meta %>%
  group_by(treatment) %>%
  summarise(n_samples = n())

# Combine and calculate prevalence
library(tidyr)  # needed for replace_na
prevalence <- entero_counts %>%
  right_join(sample_sizes, by = "treatment") %>%
  mutate(
    n_present = replace_na(n_present, 0),  # if no samples detected, set 0
    prevalence = n_present / n_samples * 100
  )

prevalence

###############################################Extended Figure 2BCD#############################################
color_code <- c(
  "Healthy"        = "#4a6325ff",
  "Compensated"    = "#f0e548ff",
  "Decompensated"  = "#e69696ff",
  "rCDI"           = "#2b6e88ff"
)

############All Taxa

Plasmid_Sum <- Type_Abundance_collapsed %>% 
  group_by(ID,treatment) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE)) %>% 
  ungroup()

Plasmid_Sum <- Plasmid_Sum %>%
  mutate(log1p_Total_TPM = log1p(Total_TPM))

Plasmid_Sum <- Plasmid_Sum %>%
  mutate(treatment = factor(treatment, 
                            levels = c("Healthy", "Compensated",  "Decompensated", "rCDI")))

Extended_Fig2B <- ggplot(Plasmid_Sum, aes(x = treatment, y = log1p_Total_TPM)) +
  geom_boxplot(aes(fill = treatment), alpha = 1, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = color_code) +
  labs(
    x = " ",
    y = "All Taxa\nPlasmid-load\nlog1p(TPM)", 
    title = " ") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =10)) +
  theme(legend.position = "none")  +
  scale_x_discrete(labels = c(
    "Healthy" = "Healthy",
    "Compensated" = "Compensated Cirrhosis",
    "Decompensated" = "Decompensated Cirrhosis",
    "rCDI" = "rCDI"))


Plasmid_Sum %>%
  filter(treatment != "Compensated") %>%
  mutate(treatment = droplevels(treatment)) %>%
  pairwise_wilcox_test(
    Total_TPM ~ treatment,
    p.adjust.method = "BH",
  ) %>%
  dplyr::select(-.y., -statistic)

############All Bacteroidaceae

Bacteroidaceae_plasmids <- Type_Abundance_collapsed %>% 
  filter(taxa_grouped =="Bacteroidaceae") %>% 
  group_by(ID,treatment) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE)) %>% 
  ungroup()

Bacteroidaceae_plasmids <- Bacteroidaceae_plasmids %>%
  mutate(log1p_Total_TPM = log1p(Total_TPM))


Bacteroidaceae_plasmids_full <- sample_meta %>%
  left_join(Bacteroidaceae_plasmids, by = c("ID","treatment")) %>%
  mutate(
    Total_TPM = replace_na(Total_TPM, 0),
    log1p_Total_TPM = replace_na(log1p_Total_TPM, 0))

Bacteroidaceae_plasmids_full <- Bacteroidaceae_plasmids_full %>%
  mutate(treatment = factor(treatment, 
                            levels = c("Healthy", "Compensated",  "Decompensated", "rCDI")))


Extended_Fig2C <-ggplot(Bacteroidaceae_plasmids_full, aes(x = treatment, y = log1p_Total_TPM)) +
  geom_boxplot(aes(fill = treatment), alpha = 1, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = color_code) +
  labs(
    x = " ",
    y = "Bacteroidaceae\nPlasmid-load\nlog1p(TPM)", 
    title = " ") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =10)) +
  theme(legend.position = "none")  +
  scale_x_discrete(labels = c(
    "Healthy" = "Healthy",
    "Compensated" = "Compensated Cirrhosis",
    "Decompensated" = "Decompensated Cirrhosis",
    "rCDI" = "rCDI"))


Bacteroidaceae_plasmids_full %>%
  filter(treatment != "Compensated") %>%
  mutate(treatment = droplevels(treatment)) %>%
  pairwise_wilcox_test(
    Total_TPM ~ treatment,
    p.adjust.method = "BH"
  ) %>%
  dplyr::select(-.y., -statistic)


############Enterobacteriaceae
Entero_plasmids <- Type_Abundance_collapsed %>% 
  filter(taxa_grouped =="Enterobacteriaceae") %>% 
  group_by(ID,treatment) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE)) %>% 
  ungroup()




Entero_plasmids <- Entero_plasmids %>%
  mutate(log1p_Total_TPM = log1p(Total_TPM))

Entero_plasmids_full <- sample_meta %>%
  left_join(Entero_plasmids, by = c("ID","treatment")) %>%
  mutate(
    Total_TPM = replace_na(Total_TPM, 0),
    log1p_Total_TPM = replace_na(log1p_Total_TPM, 0))

Entero_plasmids_full


Entero_plasmids_full <- Entero_plasmids_full %>%
  mutate(treatment = factor(treatment, 
                            levels = c("Healthy", "Compensated",  "Decompensated", "rCDI")))


Extended_Fig2D <- ggplot(Entero_plasmids_full, aes(x = treatment, y = log1p_Total_TPM)) +
  geom_boxplot(aes(fill = treatment), alpha = 1, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = color_code) +
  labs(
    x = " ",
    y = "Enterobacteriaceae\nPlasmid-load\nlog1p(TPM)", 
    title = " ") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =10)) +
  theme(legend.position = "none")  +
  scale_x_discrete(labels = c(
    "Healthy" = "Healthy",
    "Compensated" = "Compensated Cirrhosis",
    "Decompensated" = "Decompensated Cirrhosis",
    "rCDI" = "rCDI"))


Entero_plasmids_full %>%
  filter(treatment != "Compensated") %>%
  mutate(treatment = droplevels(treatment)) %>%
  pairwise_wilcox_test(
    Total_TPM ~ treatment,
    p.adjust.method = "BH"
  ) %>%
  dplyr::select(-.y., -statistic)






###### Figure

combined_plot <- Extended_Fig2A / (Extended_Fig2B  + Extended_Fig2C + Extended_Fig2D) +
  plot_layout(
    heights = c(1.5, 1)  # top panel twice the height of bottom row
  )

combined_plot
