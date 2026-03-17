###Author: Andrew J Sommer
###Date Modified: 3/16/2026
###Description: Code for Taxonomic Barplots and ANCOMBC Differential Abundance Analysis 

###Load in Library
library(ANCOMBC)
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(cowplot)
library(ggplot2)
library(patchwork)


#####Read in and Merge Metadata
sample_meta <- read.csv("Import/Metdata/Sample-Metadata.csv")
head(sample_meta)



##############################################Extended Figure 1A#############################################
Type_Abundance <- read.csv("Import/Data/Taxonomic_Abundance/Order_Microbe_Abundance_70percentcomplete_5percentcontam.csv") 
Type_Abundance <- Type_Abundance %>% left_join(sample_meta, by = "alias") #add metadata
head(Type_Abundance)


#####Deterime Relative abundance by TPM and Read Count
Type_Abundance <-
  Type_Abundance %>%
  group_by(ID) %>%                          # group by sample
  mutate(Count_Relative_Abundance = Read_Count / sum(Read_Count)) %>%  # divide by total per sample
  mutate(TPM_Relative_Abundance = TPM / sum(TPM)) %>%  # divide by total per sample
  ungroup()

Type_mean_abundance <- Type_Abundance %>%
  group_by(Taxa) %>%
  summarise(mean_Count_Relative_Abundance = mean(Count_Relative_Abundance, na.rm = TRUE),
            mean_TPM_Relative_Abundance = mean(TPM_Relative_Abundance, na.rm = TRUE)) %>%
  ungroup()

#####Deterime Low Abundance Taxa
threshold <- 0.08
low_taxa <- Type_mean_abundance  %>%
  filter(mean_Count_Relative_Abundance < threshold) %>%
  pull(Taxa)
low_taxa


#######Filter out low abundance taxa
Type_Abundance <- Type_Abundance %>%
  mutate(taxa_grouped = ifelse(Taxa %in% low_taxa, "Other", Taxa))

Type_Abundance_lumped <- Type_Abundance %>%
  # Lump all “Other” together per sample
  group_by(ID, treatment, taxa_grouped)%>%
  summarise(Count_Relative_Abundance = sum(Count_Relative_Abundance, na.rm = TRUE),
            TPM_Relative_Abundance = sum(TPM_Relative_Abundance, na.rm = TRUE))  %>%
  ungroup()




#####Set Colors and order
Taxa_colors <- c(
  "Enterobacterales" = "#CA483A",
  "Pseudomonadales"       = "#E3B44B",
  
  "Lachnospirales"    = "#7FB3E6", 
  "Lactobacillales"   = "#4A90E2",
  "Oscillospirales"   = "#1C60B3",  
  "Veillonellales"    = "#00B7B3",
  
  "Verrucomicrobiales"       = "#B07CC6", 
  
  "Bacteroidales"         = "#7B9671",
  
  "Other"               = "gray60"
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

extended_fig1A <- ggplot(Type_Abundance_lumped,
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
    axis.line = element_line(size = 0.3)) +
  theme(legend.position = "bottom")

extended_fig1A
















##############################################Extended Figure 1BC#############################################


##### Read and Clean Abundance Table. Merge with Sample data
#Abundance is from Proximeta, set to 70% completeness and 5% contamination
Abundance_Full <- read.csv("Import/Data/Taxonomic_Abundance/Family_Microbe_Abundance_70percentcomplete_5percentcontam.csv") 
Abundance_Full <- Abundance_Full %>% left_join(sample_meta, by = "alias") #add metadata
head(Abundance_Full)

##### Remove compensated group.
Abundance <- Abundance_Full %>%  filter(treatment != "Compensated")



##### Create a phyloseq object based on read count
# ANCOM-BC works on count data and has its own normalization procedure

otu_table <- Abundance %>%
  dplyr::select(ID, Family, Read_Count) %>%
  tidyr::pivot_wider(names_from = ID, values_from = Read_Count, values_fill = 0) %>%
  column_to_rownames("Family") %>%
  as.matrix()
otu_table <- otu_table(otu_table, taxa_are_rows = TRUE)

sample_data <- Abundance %>%
  dplyr::select(ID, treatment) %>%
  distinct() %>%
  column_to_rownames("ID") %>%
  as.data.frame()
sample_data <- sample_data(sample_data)

ps <- phyloseq(otu_table, sample_data)
ps



##### Prepare data for ANCOM-BC
#Remove Low abundance families
ps_filtered <- filter_taxa(ps, function(x) sum(x > 0) > 0.10 * nsamples(ps), TRUE)
ntaxa(ps)
ntaxa(ps_filtered) 

#Re-level so that healthy is baseline
sample_data(ps_filtered)$treatment <- relevel(factor(sample_data(ps_filtered)$treatment),ref = "Healthy")
levels(sample_data(ps_filtered)$treatment)


##### Run ANCOM-BC
out_ancombc <- ancombc(
  data = ps_filtered,
  formula = "treatment",    
  p_adj_method = "BH",
  lib_cut = 0,
  group = "treatment",
  verbose = TRUE
)


##### Extract ANCOM-BC results
res <- out_ancombc$res

#Extract and clean the log fold change (effect size) data
da_lfc <- res$lfc %>% dplyr::select(-`(Intercept)`) %>%
  rename(
    lfc_Decompensated  = treatmentDecompensated,
    lfc_rCDI = treatmentrCDI
  )

#Extract and clean the adjusted p values (BH/FDR) 
da_pval <- res$q_val %>% dplyr::select(-`(Intercept)`) %>%
  rename(
    BH_Decompensated  = treatmentDecompensated,
    BH_rCDI = treatmentrCDI
  )

ancombc_output <- da_lfc %>%
  left_join(da_pval, by = "taxon")

head(ancombc_output)
unique(ancombc_output$taxon)


##### Decompensated Group Signficant Taxa
sig_Decomp <- ancombc_output %>%
  filter(BH_Decompensated < 0.05) %>%
  dplyr::select(taxon, lfc_Decompensated, BH_Decompensated)

sig_Decomp <- sig_Decomp %>%
  mutate(
    taxon = case_when(
      taxon == "CAG-508" ~ "CAG-508 (Clostridia)",
      taxon == "UBA1381" ~ "UBA1381 (Clostridia)",
      TRUE ~ taxon))

sig_Decomp <- sig_Decomp %>%
  arrange(lfc_Decompensated) %>%
  mutate(taxon = factor(taxon, levels = taxon))


sig_Decomp

##### rCDI Group Signficant Taxa
sig_rCDI <- ancombc_output %>%
  filter(BH_rCDI < 0.05) %>%
  dplyr::select(taxon, lfc_rCDI, BH_rCDI)

sig_rCDI <- sig_rCDI %>%
  mutate(
    taxon = case_when(
      taxon == "CAG-508" ~ "CAG-508 (Clostridia)",
      taxon == "ICN-92133" ~ "ICN-92133 (Negativicutes)",
      taxon == "CAG-508" ~ "CAG-508 (Clostridia)",
      TRUE ~ taxon))

sig_rCDI <- sig_rCDI %>%
  arrange(lfc_rCDI) %>%
  mutate(taxon = factor(taxon, levels = taxon))


sig_rCDI


##### Plot  Signficant Taxa

extended_fig1B <- ggplot(sig_Decomp, aes(x = taxon, y = lfc_Decompensated, fill = lfc_Decompensated > 0)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +  
  coord_flip() +  # horizontal bars
  scale_fill_manual(
    values = c("#4a6325ff", "#e69696ff"),
    name = "Association",
    labels = c("Healthy", "Decompensated Cirrhosis")
  ) +
  labs(x = " ", y = "Log Fold Change (Effect Size)", title = " ") +
  theme_classic(base_size = 11) +  # base font size 11
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    plot.title = element_text(size = 11)
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)
extended_fig1B

extended_fig1C <- ggplot(sig_rCDI, aes(x = taxon, y = lfc_rCDI, fill = lfc_rCDI > 0)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +  
  coord_flip() +  # horizontal bars
  scale_fill_manual(
    values = c("#4a6325ff", "#2b6e88ff"),
    name = "Association",
    labels = c("Healthy", "rCDI")
  ) +
  labs(x = " ", y = "Log Fold Change (Effect Size)", title = " ") +
  theme_classic(base_size = 11) +  # set base font size
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)
extended_fig1C

