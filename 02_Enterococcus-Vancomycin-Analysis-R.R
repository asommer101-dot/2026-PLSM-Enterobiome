###Author: Andrew J Sommer
###Date Modified: 3/17/2026
###Description: Code for Anaylsis of Enterococcus and Vancomycin ARGs
library(dplyr)
library(vegan)
library(tibble) 
library(ggplot2)
library(MetBrewer)
library(cowplot)
library(forcats)
library(stringr)
library(tidyr)
library(rstatix)

#####Load in metadata
sample_meta <- read.csv("Import/Metadata/Sample-Metadata.csv")
head(sample_meta)


###############################################FIGURE 2B#############################################

#####Load in ARG Table 
arg_cleaned <- read.csv("Import/Data/AMR_VF/arg-table-core.csv")

#####Subset for Enterococcus Vancomycin (Glycopeptide) Genes
Enterococcus_ARGs <- arg_cleaned %>%  filter(Family == "Enterococcaceae")
Enterococcus_vanco <- Enterococcus_ARGs%>%  filter(Gene_class == "GLYCOPEPTIDE")
head(Enterococcus_vanco)

#####Summarize MAG Enterococcus Vancomycin carriage
genome_gene_summary <- Enterococcus_vanco %>%
  group_by(ID, Genome_ID, gene_id, Contig_type) %>%
  summarise(
    Contig_TPM = mean(Contig_TPM, na.rm = TRUE),
    Genome_TPM = mean(Genome_TPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(log1p_Contig_TPM = log1p(Contig_TPM)) %>%
  mutate(log1p_Genome_TPM = log1p(Genome_TPM)) 


genome_gene_summary <- genome_gene_summary %>%
  mutate(
    Genome_Simple = Genome_ID %>%
      gsub("CD04_bin_6", "CD04 (Bin_6)", .) %>%
      gsub("CD13_bin_12", "CD13 (Bin_12)", .) %>%
      gsub("CD13_bin_16", "CD13 (Bin_16)", .) %>%
      gsub("CD17_bin_13", "CD17 (Bin_13)", .) %>%
      gsub("CD17_bin_46", "CD17 (Bin_46)", .) %>%
      gsub("CD29_bin_7", "CD29 (Bin_7)", .) %>%
      gsub("CD32_bin_45", "CD32 (Bin_45)", .) %>%
      gsub("LC14_bin_21", "LC14 (Bin_21)", .) %>%
      gsub("LC19_bin_12", "LC19 (Bin_12)", .))

#####Plot Fig2B
Fig2B <- ggplot(genome_gene_summary, aes(x = gene_id, y = Genome_Simple, fill = log1p_Contig_TPM)) +
  geom_tile(color = "black") +  # black borders
  scale_fill_gradient(low = "white", high = "peru") +
  theme_classic(base_size = 11) +
  geom_text(
    data = genome_gene_summary %>% filter(Contig_type == "plasmid"),
    aes(label = "P"),
    color = "black",
    fontface = "bold",
    size = 3
  )  +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.position = "top"
  ) +
  labs(
    x = "",
    y = " ",
    fill = "log1p(Contig TPM)") 
Fig2B






###############################################Extended Fig5#############################################


#####Load in ARG Table 
arg_cleaned <- read.csv("Import/Data/AMR_VF/arg-table-core.csv")

#####Subset to Vancomycin ARGs 
Vancomycin_ARGs <- arg_cleaned %>%  filter(Gene_class == "GLYCOPEPTIDE")
Vancomycin_ARGs <- Vancomycin_ARGs %>% 
  dplyr::select(ID,treatment,Class,Order,Family,Genus, Contig_type, Contig_TPM)

head(Vancomycin_ARGs)


############################# Genomic ARGs (Vancomycin)
Vancomycin_gARGs_TPM <- Vancomycin_ARGs %>% 
  filter(Contig_type == "genomic")
head(Vancomycin_gARGs_TPM)
unique(Vancomycin_gARGs_TPM$ID)
sum(Vancomycin_gARGs_TPM$Contig_TPM)


Vancomycin_gARGs_TPM %>% 
  group_by(Class,Order,Family,Genus,Contig_type)  %>%  
  summarise(
    total_TPM = sum(Contig_TPM, na.rm = TRUE),
    .groups = "drop"
  )

Vancomycin_gARGs_TPM <- Vancomycin_gARGs_TPM %>%
  mutate(
    TaxGroup = case_when(
      Genus %in% c("Enterococcus_D") ~ "Enterococcus",
      Genus %in% c("Paenibacillus_B") ~ "Paenibacillaceae",
      Genus %in% c("Clostridioides") ~ "Clostridioides",
      Genus %in% c("Butyricicoccus") ~ "Butyricicoccaceae",
      Genus %in% c("Clostridium_J") ~ "Clostridiaceae",
      Genus %in% c("Hominenteromicrobium", "Clostridium_A") ~ "Acutalibacteraceae",
      Genus %in% c("Enterocloster","Cuneatibacter","Luxibacter","Bariatricus","Clostridium_AP","Blautia","Anaerostipes","Blautia_A","Fimimorpha") ~ "Lachnospiraceae",
      Genus %in% c("Massilioclostridium", "Massiliimalia") ~ "Ruminococcaceae"))


Vancomycin_gARGs_TPM_Summary <-  Vancomycin_gARGs_TPM %>%
  group_by(ID, TaxGroup, treatment, Contig_type) %>%   # include treatment if you want to keep facets
  summarise(
    total_TPM = sum(Contig_TPM, na.rm = TRUE),
    .groups = "drop")

tax_colors <- c(
  "Acutalibacteraceae"    = "#66c2a5", 
  "Butyricicoccaceae"     = "#41b6c4",  
  "Clostridiaceae"        = "#1f78b4",  
  "Lachnospiraceae"       = "#a6cee3",  
  "Ruminococcaceae"       = "#b2df8a",  
  "Paenibacillaceae" = "#9EFAB8", 
  "Enterococcus"       = "#ff7f00",  
  "Clostridioides"      = "#fdbf6f")


Vancomycin_gARGs_TPM_Summary <- Vancomycin_gARGs_TPM_Summary %>%
  mutate(TaxGroup = factor(TaxGroup, 
                           levels = c("Clostridioides","Enterococcus",
                                      "Acutalibacteraceae","Butyricicoccaceae","Clostridiaceae",
                                      "Lachnospiraceae","Ruminococcaceae","Paenibacillaceae")))   



Vancomycin_gARGs_TPM_Summary <- Vancomycin_gARGs_TPM_Summary %>%
  mutate(log1p_TPM = log1p(total_TPM))



Vancomycin_gARGs_TPM_Summary$ID <- factor(
  Vancomycin_gARGs_TPM_Summary$ID,
  levels = c(
    sort(grep("^LC", unique(Vancomycin_gARGs_TPM_Summary$ID), value = TRUE)),
    sort(grep("^CD", unique(Vancomycin_gARGs_TPM_Summary$ID), value = TRUE))
  )
)


gARG_VRE <- ggplot(Vancomycin_gARGs_TPM_Summary, aes(x = ID, y = log1p_TPM, fill = TaxGroup)) +
  geom_bar(stat = "identity", color = "black", width = 1, linewidth = 0.05) +  
  theme_half_open() +
  scale_y_continuous( limits = c(0, 20), expand = c(0, 0) ) +
  labs(
    x = "Sample",
    y = "Chromosomal VRGs\n log1p(TPM)",
    fill = "Family",
  ) +
  scale_fill_manual(values = tax_colors) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  ) +
  theme(legend.position = "none")
gARG_VRE


Vancomycin_ARGs_ContigTypes_summary <- Vancomycin_ARGs %>%
  group_by(Contig_type) %>%
  summarise(
    total_TPM = sum(Contig_TPM, na.rm = TRUE),   # or Contig_TPM if you prefer
    count = n(),                                 # number of rows per contig type
    .groups = "drop"
  ) %>%
  mutate(log1p_TPM = log1p(total_TPM)) 





############################# Plasmid ARGs (Vancomycin)
Vancomycin_pARGs_TPM <- Vancomycin_ARGs %>% 
  filter(Contig_type == "plasmid")
head(Vancomycin_pARGs_TPM)
unique(Vancomycin_pARGs_TPM$ID)
sum(Vancomycin_pARGs_TPM$Contig_TPM)
unique(Vancomycin_pARGs_TPM$Genus)


Vancomycin_pARGs_TPM <- Vancomycin_pARGs_TPM %>%
  mutate(
    TaxGroup = case_when(
      Genus %in% c("Enterococcus_D") ~ "Enterococcus",
      Genus %in% c("Unclassified_Bacteria") ~ "Unclassified"))


unique(Vancomycin_pARGs_TPM$TaxGroup)

Vancomycin_pARGs_TPM_Summary <-  Vancomycin_pARGs_TPM %>%
  group_by(ID, TaxGroup, treatment, Contig_type) %>%   # include treatment if you want to keep facets
  summarise(
    total_TPM = sum(Contig_TPM, na.rm = TRUE),
    .groups = "drop")

tax_colors_2 <- c(
  "Enterococcus"       = "#ff7f00",  
  "Unclassified"      = "grey80")

Vancomycin_pARGs_TPM_Summary <- Vancomycin_pARGs_TPM_Summary %>%
  mutate(log1p_TPM = log1p(total_TPM))

pARG_VRE <- ggplot(Vancomycin_pARGs_TPM_Summary, aes(x = ID, y = log1p_TPM, fill = TaxGroup)) +
  geom_bar(stat = "identity", color = "black", width = 1, linewidth = 0.05) +  
  theme_half_open() +
  scale_y_continuous( limits = c(0, 20), expand = c(0, 0) ) +
  labs(
    x = "Sample",
    y = "Plasmid VRGs\n log1p(TPM)",
    fill = "Family",
  ) +
  scale_fill_manual(values = tax_colors_2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")

VRE_Combined <- plot_grid(pARG_VRE, gARG_VRE,rel_widths = c(1,4))
VRE_Combined



