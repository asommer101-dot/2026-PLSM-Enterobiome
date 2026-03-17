###Author: Andrew J Sommer
###Date Modified: 3/16/2026
###Description: Code for Anaylsis of Enterobactericeae ARGs and VFs
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

#####Set up color codes
color_code <- c(
  "Healthy"        = "#4a6325ff",
  "Compensated"    = "#f0e548ff",
  "Decompensated"  = "#e69696ff",
  "rCDI"           = "#2b6e88ff")

family_colors <- c("Enterobacteriaceae" = "#CA483A",
                   "Other" = "grey60",
                   "Unclassified" = "grey80")


#####Load in metadata
sample_meta <- read.csv("Import/Metadata/Sample-Metadata.csv")
head(sample_meta)



###############################################FIGURE 1#############################################
#####Load in ARG Table
arg_cleaned <- read.csv("Import/Data/AMR_VF/arg-table-core.csv")

#####Fix ARG Names to make plotting more intuitive
#The PHENICOL/QUINOLONE are oqx efflux pumps
#Combine MACROLIDE/LINCOSAMIDE related classes together
arg_cleaned <- arg_cleaned %>%
  mutate(Gene_class = recode(Gene_class,
                             "AMINOGLYCOSIDE/QUINOLONE" = "AMINOGLYCOSIDE",
                             "LINCOSAMIDE" = "MACROLIDE/LINCOSAMIDE",
                             "MACROLIDE" = "MACROLIDE/LINCOSAMIDE",
                             "MACROLIDE/LINCOSAMIDE/STREPTOGRAMIN" = "MACROLIDE/LINCOSAMIDE",
                             "LINCOSAMIDE/STREPTOGRAMIN" = "MACROLIDE/LINCOSAMIDE",
                             "PHENICOL/QUINOLONE" = "EFFLUX")) 
#write.csv(arg_cleaned, "Import/Data/AMR_VF/arg-table-cleaned.csv")



#####Remove low abudance ARG Classes
#Summarize Unique Hits in each class
arg_cleaned %>%
  group_by(Gene_class) %>%
  summarise(
    total_hits = n(),
    unique_genes = n_distinct(gene_id)
  ) %>%
  arrange(desc(unique_genes))

#Remove the following gene classes given lowest overall abundance which would clutter plots
arg_filtered <- arg_cleaned %>%
  filter(!Gene_class %in% c("PLEUROMUTILIN", "STREPTOTHRICIN", "COLISTIN",
                            "BLEOMYCIN", "FUSIDIC ACID", "MUPIROCIN", "RIFAMYCIN",
                            "SULFONAMIDE"))
arg_filtered <- arg_filtered %>% mutate(Gene_class = str_to_title(Gene_class))

#Calculate Gene counts in the filtered set
gene_counts <- arg_filtered %>%
  group_by(Gene_class) %>%
  summarise(
    total_hits = n(),
    unique_genes = n_distinct(gene_id)) 

#####Fix naming schme for Unclassified Genomes
arg_filtered <- arg_filtered %>%
  mutate(Family = replace(Family, Family == "Unclassified_Bacteria", "Unclassified"),
         Family = replace(Family, is.na(Family), "Unclassified"))

arg_filtered <- arg_filtered %>%
  mutate(Order = replace(Order, Order == "Unclassified_Bacteria", "Unclassified"),
         Order = replace(Order, is.na(Order), "Unclassified"))

arg_filtered <- arg_filtered %>%
  mutate(Class = replace(Class, Class == "Unclassified_Bacteria", "Unclassified"),
         Class = replace(Class, is.na(Class), "Unclassified"))



#####Determine Relative Abundance of Enterobactericeae ARGs based off TPM (Fig1A)

#Relative contribution of bacterial families encoding each gene class (TPM)
arg_stack <- arg_filtered %>%
  group_by(Gene_class, Family,Order,Class) %>%
  summarise(TPM = sum(Contig_TPM, na.rm = TRUE), .groups = "drop") %>%
  group_by(Gene_class) %>%
  mutate(rel_abundance = TPM / sum(TPM)) %>%
  ungroup()

#Group as either Enterobacteriaceae, Unclassified, or Other
arg_stack_grouped <- arg_stack %>%
  mutate(
    Family_group = case_when(
      Family == "Enterobacteriaceae" ~ "Enterobacteriaceae",
      Family == "Unclassified" ~ "Unclassified",
      TRUE ~ "Other")) %>%
  group_by(Gene_class, Family_group) %>%
  summarise(
    sum_rel_abundance = sum(rel_abundance, na.rm = TRUE),.groups = "drop") 

#Merge with the counts of genes
arg_stack_grouped <- arg_stack_grouped %>% left_join(gene_counts, by = "Gene_class")


#Plot the relative contribution (TPM) of ARGs by Enterobacteriaceae
arg_stack_grouped <- arg_stack_grouped %>%
  mutate(Gene_class = factor(Gene_class, levels = sort(unique(Gene_class), decreasing = TRUE)))  # reverse alphabetical

Fig1A <- ggplot(arg_stack_grouped, 
  aes(x = Gene_class, y = sum_rel_abundance, fill = Family_group)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = family_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  
  coord_flip() +  # make bars horizontal
  theme_minimal(base_size = 10) +
  theme_half_open() +
  theme(
    axis.text.y = element_text(size = 10),       
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1.2, "lines")) +
  labs(
    x = "",
    y = "Relative Abundance",
    fill = " ") +
  geom_text(
    data = distinct(arg_stack_grouped, Gene_class, unique_genes),
    aes(x = Gene_class, y = 1.05, label = paste0("n=", unique_genes)),
    size = 3,
    fontface = "bold",
    inherit.aes = FALSE) +
  guides(fill = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom")


#####Subset ARGs carried by Enterobacteriaceae on chromosomes (Fig1B) or plasmids (Fig1C)
ARG_entero <- arg_cleaned %>%  filter(Family == "Enterobacteriaceae") #arg_cleaned contains the low diversity ARG groups, which we retain for this analysis
ARG_entero_plasmid <- ARG_entero %>%  filter(Contig_type == "plasmid")
ARG_entero_genomic <- ARG_entero %>%  filter(Contig_type == "genomic")



#####Chromosomal-Borne Enterobacteriaceae ARGs

summary_entero_genomic <- ARG_entero_genomic %>%
  group_by(ID) %>%
  summarise(Sum_Contig_TPM = sum(Contig_TPM, na.rm = TRUE)) %>%
  complete(ID = unique(sample_meta$ID), fill = list(Sum_Contig_TPM = 0)) %>%
  mutate(log1p_Sum = log1p(Sum_Contig_TPM))

summary_entero_genomic <- sample_meta %>%
  left_join(summary_entero_genomic, by = "ID")

summary_entero_genomic <- summary_entero_genomic %>%
  mutate(treatment = factor(treatment, 
                            levels = c("Healthy", "Compensated",  "Decompensated", "rCDI")))

#Pairwise Wilcox to compare Chromosomal-borne ARG abundance to rCDI 
summary_entero_genomic %>%
  filter(treatment != "Compensated") %>%
  mutate(treatment = droplevels(treatment)) %>%
  pairwise_wilcox_test(
    Sum_Contig_TPM ~ treatment,
    p.adjust.method = "BH",
  ) %>%
  dplyr::select(-.y., -statistic, -p)

#Pairwise Wilcox to compare Chromosomal-borne ARG abundance to Healthy 
summary_entero_genomic %>%
  filter(treatment != "Compensated") %>%
  mutate(treatment = droplevels(treatment)) %>%
  pairwise_wilcox_test(
    Sum_Contig_TPM ~ treatment,
    p.adjust.method = "BH",
  ) %>%
  dplyr::select(-.y., -statistic, -p)

Fig1B <- ggplot(summary_entero_genomic, aes(x = treatment, y = log1p_Sum)) +
  geom_boxplot(aes(fill = treatment), alpha = 1, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = color_code) +
  labs(
    x = " ",
    y = "Enterobacteriaceae\nChromosomal ARGs\nlog1p(TPM)", 
    title = " ") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =10)) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Healthy" = "Healthy",
    "Compensated" = "Compensated Cirrhosis",
    "Decompensated" = "Decompensated Cirrhosis",
    "rCDI" = "rCDI"))




#####Plasmid-Borne Enterobacteriaceae ARGs
summary_entero_plasmid <- ARG_entero_plasmid %>%
  group_by(ID) %>%
  summarise(Sum_Contig_TPM = sum(Contig_TPM, na.rm = TRUE)) %>%
  complete(ID = unique(sample_meta$ID), fill = list(Sum_Contig_TPM = 0)) %>%
  mutate(log1p_Sum = log1p(Sum_Contig_TPM))

summary_entero_plasmid <- sample_meta %>%
  left_join(summary_entero_plasmid, by = "ID")

summary_entero_plasmid <- summary_entero_plasmid %>%
  mutate(treatment = factor(treatment, 
                            levels = c("Healthy", "Compensated",  "Decompensated", "rCDI")))

#Pairwise Wilcox to compare plasmid-borne ARG abundance to rCDI 
summary_entero_plasmid %>%
  pairwise_wilcox_test(
    Sum_Contig_TPM ~ treatment,
    p.adjust.method = "BH",
    ref.group = "rCDI") %>%
  dplyr::select(-.y., -statistic, -p)

#Pairwise Wilcox to compare plasmid-borne ARG abundance to Healthy
summary_entero_plasmid %>%
  pairwise_wilcox_test(
    Sum_Contig_TPM ~ treatment,
    p.adjust.method = "BH",
    ref.group = "Healthy") %>%
  dplyr::select(-.y., -statistic, -p)

#Plot the boxplot
Fig1C <- ggplot(summary_entero_plasmid, aes(x = treatment, y = log1p_Sum)) +
  geom_boxplot(aes(fill = treatment), alpha = 1, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = color_code) +
  labs(
    x = " ",
    y = "Enterobacteriaceae\nPlasmid-borne ARGs\nlog1p(TPM)", 
    title = " ") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =10)) +
  theme(legend.position = "none")  +
  scale_x_discrete(labels = c(
    "Healthy" = "Healthy Control",
    "Compensated" = "Compensated Cirrhosis",
    "Decompensated" = "Decompensated Cirrhosis",
    "rCDI" = "rCDI"))



#####Summarize Enterobacteriaceae Virulence Factor Load 
vf_table_cleaned <- read.csv(file = "Import/Data/AMR_VF/virulence-gene-table.csv") #cov.adj > 90%
vf_table_cleaned %>% count(Family) #Nearly all VFs (n=379) originate from Enterobacteriaceae. The remaining unclassified VFs may also be Enterobacteriaceae, but are not linked to a genome (even a low-confidence genome) and therefore excluded from this analysis


vf_entero <- vf_table_cleaned %>%  filter(Family == "Enterobacteriaceae")

summary_entero_vf <- vf_entero %>%
  complete(ID = sample_meta$ID) %>%  
  group_by(ID) %>%
  summarise(VF_TPM = sum(Contig_TPM, na.rm = TRUE)) %>%
  mutate(VF_log1p = log1p(VF_TPM)) %>%
  left_join(sample_meta, by = "ID")


#####Scatterplot of Enterobacteriaceae Virulence Factor and ARG Abundance
summary_entero_arg <- ARG_entero %>%
  complete(ID = sample_meta$ID) %>%  
  group_by(ID) %>%
  summarise(ARG_TPM = sum(Contig_TPM, na.rm = TRUE)) %>%
  mutate(ARG_log1p = log1p(ARG_TPM)) %>%
  left_join(sample_meta, by = "ID") %>%
  dplyr::select(-treatment, -alias)

scatter_df <- summary_entero_arg %>%
  left_join(summary_entero_vf, by = c("ID"))

scatter_df <- scatter_df %>%
  mutate(treatment = factor(treatment, 
                            levels = c("Healthy", "Compensated",  "Decompensated", "rCDI")))


#Plot ARG vs VF Abundance
Fig1D <- ggplot(scatter_df, aes(x = ARG_log1p, y = VF_log1p)) +
  theme_classic() +
  geom_point(aes(fill = treatment),  
             shape = 21,             
             color = "black",        
             size = 3.5,               
             stroke = 0.4) +        
  scale_fill_manual(values = color_code) +  
  theme(legend.position = "none") +
  labs(
    x = "Enterobacteriaceae\nARGs log1p(TPM)",
    y = "Enterobacteriaceae\nVFs log1p(TPM)",
    title = "",
    fill = "Treatment")

#Plot ARG vs VF Abundance. Jittered to see (0,0) better
Fig1D_Jitter <- ggplot(scatter_df, aes(x = ARG_log1p, y = VF_log1p)) +
  theme_classic() +
  geom_jitter(
    aes(fill = treatment),
    shape = 21,
    color = "black",
    size = 3.5,
    stroke = 0.4,
    width = 0.1,
    height = 0.1) +        
  scale_fill_manual(values = color_code) +  
  theme(legend.position = "none") +
  labs(
    x = "Enterobacteriaceae\nARGs log1p(TPM)",
    y = "Enterobacteriaceae\nVFs log1p(TPM)",
    title = "",
    fill = "Treatment")

cor.test(scatter_df$ARG_log1p, scatter_df$VF_log1p, method = "spearman")


#####Combine to create Figure 1 Panels. Note that figure is realigned in graphics editor
left_short <- Fig1A
right_stack <- plot_grid(Fig1B, Fig1C, Fig1D, ncol = 1, align = "v")

combined_plot <- plot_grid(
  left_short, right_stack,
  ncol = 2,
  rel_widths = c(1.5, 1))
combined_plot



############################################### Extended_Fig4 (MAG Level Scatterplot) #############################################

#####Calculate ARG counts within each Genome bin
Entero_ARG_counts <- ARG_entero %>%
  group_by(Genome_ID, Genus, Species, treatment, ID, Contig_type) %>%
  summarise(count = n_distinct(gene_id), .groups = "drop") %>%
  pivot_wider(
    names_from = Contig_type,
    values_from = count,
    values_fill = 0
  ) %>%
  rename(ARG_Genomic = genomic,
         ARG_Plasmid = plasmid) %>%
  mutate(ARG_Total = ARG_Genomic + ARG_Plasmid)

#####Calculate VF counts within each Genome bin
Entero_VF_counts <- vf_entero %>%
  group_by(Genome_ID, Genus, Species, treatment, ID, Contig_type) %>%
  summarise(count = n_distinct(gene_id), .groups = "drop") %>%
  pivot_wider(
    names_from = Contig_type,
    values_from = count,
    values_fill = 0
  ) %>%
  rename(VF_Genomic = genomic,
         VF_Plasmid = plasmid,
         VF_Viral = viral) %>%
  mutate(VF_Total = VF_Genomic + VF_Plasmid + VF_Viral)


Entero_genomes_merged <- full_join(
  Entero_ARG_counts,
  Entero_VF_counts,
  by = c("Genome_ID","Genus", "Species", "treatment", "ID")) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))
  
#####Merge in Enterobacteriaceae MAG Data
#Enterobacteriaceae-MAGs-Full.csv was generated by sub-setting the full MAG set, which is a large file. It contains all MAGs regardless of completeness or contamination
#This step is optional, but helpful if you want to further subset plots by completeness or contamination levels
MAGs_Entero <- read.csv("Import/Data/MAGs/Enterobacteriaceae-MAGs-Full.csv")
MAGs_Entero <- MAGs_Entero %>% 
  dplyr::select(genome_bin,Contamination,Completeness) %>% 
  dplyr::rename(Genome_ID = genome_bin)

Entero_genomes_merged_complete <- Entero_genomes_merged %>%
  left_join(MAGs_Entero, by = "Genome_ID")


#####Subset Escherichia and Plot Extended_Fig4A
Escherichia_Genomes <- Entero_genomes_merged_complete %>%
  filter(Genus == "Escherichia") %>% 
  dplyr::select(Genome_ID, ID, Species, treatment, ARG_Total, VF_Total, Contamination, Completeness)

Escherichia_Genomes <- Escherichia_Genomes %>%
  mutate(treatment = factor(treatment, 
                            levels = c("Healthy", "Compensated",  "Decompensated", "rCDI")))

cor.test(Escherichia_Genomes$ARG_Total, Escherichia_Genomes$VF_Total, method = "spearman")

Extended_Fig4A <- ggplot(Escherichia_Genomes, aes(x = ARG_Total, y = VF_Total)) +
  theme_classic() +
  geom_point(aes(fill = treatment),  
             shape = 21,             
             color = "black",        
             size = 3.5,               
             stroke = 0.1,
             position = position_jitter(width = 0.2, height = 0.2)) +  # add jitter
  scale_fill_manual(values = color_code) +  
  theme(legend.position = "none") +
  labs(
    x = "Number ARGs",
    y = "Number VFs",
    fill = "Treatment")

Extended_Fig4A

#####Subset Klebsiella and Plot Extended_Fig4B
Klebsiella_Genomes <- Entero_genomes_merged_complete %>%
  filter(Genus == "Klebsiella") %>% 
  dplyr::select(Genome_ID, ID, Species, treatment, ARG_Total, VF_Total, Contamination, Completeness)

Klebsiella_Genomes <- Klebsiella_Genomes %>%
  mutate(treatment = factor(treatment, 
                            levels = c("Healthy", "Compensated",  "Decompensated", "rCDI")))

cor.test(Klebsiella_Genomes$ARG_Total, Klebsiella_Genomes$VF_Total, method = "spearman")


Extended_Fig4B <- ggplot(Klebsiella_Genomes, aes(x = ARG_Total, y = VF_Total)) +
  theme_classic() +
  geom_point(aes(fill = treatment),  
             shape = 21,             
             color = "black",        
             size = 3.5,               
             stroke = 0.1,
             position = position_jitter(width = 0.2, height = 0.2)) +  # add jitter
  scale_fill_manual(values = color_code) +  
  theme(legend.position = "none") +
  labs(
    x = "Number ARGs",
    y = "Number VFs",
    fill = "Treatment"
  )
Extended_Fig4B


############################################### Extended_Fig3 (Heatmap of Enterobacteriaceae ARGs) #############################################

#####Chromosomal ARG Heatmap (Not Included in Manuscript)
ARG_entero_genomic_2 <- ARG_entero_genomic %>%
  dplyr::mutate(
    treatment = factor(treatment, 
                       levels = c("Healthy", "Compensated", "Decompensated", "rCDI")),
    sample = fct_reorder(ID, as.numeric(treatment)))

ARG_entero_genomic_2 <- ARG_entero_genomic_2 %>% 
  dplyr::mutate(gene_id = factor(gene_id, levels = unique(gene_id[order(Gene_class)])))

Entero_Genomic_long <- ARG_entero_genomic_2 %>%
  group_by(ID, gene_id) %>%
  summarise(total_TPM = sum(Contig_TPM, na.rm = TRUE), .groups = "drop") %>%
  complete(ID, gene_id, fill = list(total_TPM = 0))


Entero_Genomic_long <- Entero_Genomic_long %>%
  left_join(ARG_entero_genomic_2 %>% distinct(ID, treatment), by = "ID") %>%
  left_join(ARG_entero_genomic_2 %>% distinct(gene_id, Gene_class),by = "gene_id")


###Samples WITH Entero Genomic-borne ARGs###
Entero_Genomic_long  %>%
  select(ID, treatment) %>%
  unique()


Entero_Genomic_long <- Entero_Genomic_long %>% left_join(sample_meta, by = "ID")

Entero_Genomic_long$ID <- factor(
  Entero_Genomic_long$ID,
  levels = c(
    sort(grep("^HC", unique(Entero_Genomic_long$ID), value = TRUE)),
    sort(grep("^LC", unique(Entero_Genomic_long$ID), value = TRUE)),
    sort(grep("^CD", unique(Entero_Genomic_long$ID), value = TRUE))
  )
)


Entero_Genomic_long <- Entero_Genomic_long %>%
  # convert gene_id to character for proper alphabetical sorting
  mutate(gene_id = as.character(gene_id)) %>%
  arrange(Gene_class, gene_id) %>%
  # set factor levels in the sorted order
  mutate(gene_id = factor(gene_id, levels = unique(gene_id)))


ggplot(Entero_Genomic_long, aes(x = ID, y = gene_id))  +
  geom_tile(aes(fill = total_TPM), color = "black", linewidth = 0.05) +
  geom_point(
    data = Entero_Genomic_long %>% filter(total_TPM > 0),
    aes(x = ID, y = gene_id),
    color = "black",
    size = 1
  ) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(x = " ", y = "", fill = "TPM") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=10),
    axis.text.y = element_text(size=10)  )



#####Plasmid ARG Heatmap (Extended_Fig3)
Entero_ARG_Plasmid_long <- ARG_entero_plasmid %>%
  group_by(ID, gene_id) %>%
  summarise(total_TPM = sum(Contig_TPM, na.rm = TRUE), .groups = "drop") %>%
  complete(ID, gene_id, fill = list(total_TPM = 0))


Entero_ARG_Plasmid_long <- Entero_ARG_Plasmid_long %>%
  left_join(ARG_entero_plasmid %>% distinct(ID, treatment), by = "ID") %>%
  left_join(ARG_entero_plasmid %>% distinct(gene_id, Gene_class),by = "gene_id")

Entero_Genus_Prop_ARG <- ARG_entero_plasmid %>%
  group_by(gene_id, Genus) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(
    prop = n / sum(n)   # proportion of each species within this gene_id
  ) %>%
  ungroup()

Entero_Genus_Prop_ARG <- Entero_Genus_Prop_ARG %>% mutate(Genus = recode(Genus, "Citrobacter_A" = "Citrobacter"))

genus_colors <- c(
  "Klebsiella"   = "#6b200c",
  "Escherichia"  = "#973d21",
  "Morganella"   = "#da6c42",
  "Citrobacter"  = "#ee956a",
  "Proteus"      = "#fbc2a9",
  "Kluyvera"     = "#f6f2ee"
)


Entero_ARG_Plasmid_long <- Entero_ARG_Plasmid_long %>% left_join(sample_meta, by = "ID")

Entero_ARG_Plasmid_long$ID <- factor(
  Entero_ARG_Plasmid_long$ID,
  levels = c(
    sort(grep("^HC", unique(Entero_ARG_Plasmid_long$ID), value = TRUE)),
    sort(grep("^LC", unique(Entero_ARG_Plasmid_long$ID), value = TRUE)),
    sort(grep("^CD", unique(Entero_ARG_Plasmid_long$ID), value = TRUE))
  )
)

Plasmid_ARG_PlotA <- ggplot(Entero_ARG_Plasmid_long, aes(x = ID, y = gene_id))  +
  geom_tile(aes(fill = total_TPM), color = "black", linewidth = 0.05) +
  geom_point(
    data = Entero_ARG_Plasmid_long %>% filter(total_TPM > 0),
    aes(x = ID, y = gene_id),
    color = "black",
    size = 1
  ) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(x = " ", y = "", fill = "TPM") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=10),
    axis.text.y = element_text(size=10)  )


Plasmid_ARG_PlotB <-ggplot(Entero_Genus_Prop_ARG, aes(x = gene_id, y = prop, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.1) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    legend.position = "right" ) +
  coord_flip() +
  theme(
    axis.title.y = element_blank(),  
    axis.title.x = element_blank(),    
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=10)) +
  scale_fill_manual(values = genus_colors)

library(patchwork)
Extended_Fig3 <- Plasmid_ARG_PlotA + Plasmid_ARG_PlotB + plot_layout(widths = c(2, 1)) + plot_layout(guides = "collect") & theme(legend.position = "right")
Extended_Fig3


###############################################Fig2A High Risk Enterobacteriaceae MAGs#############################################


#####Code to plot Enterobacteriaceae MAGs with plasmid borne beta-lactamases
Plasmid_Source <- ARG_entero_plasmid %>%
  select(gene_id, Genome_ID, Genome_TPM, Genus, Species, ID, Source_Genome, treatment) %>%
  distinct()  %>%
  group_by(ID, gene_id, Genus, Species, Genome_TPM, Genome_ID,Source_Genome, treatment) %>%  # add metadata columns
  summarise(Count = n(), .groups = "drop") %>% 
  rename(Present = Count)

Plasmid_Source <- Plasmid_Source %>%
  group_by(Genome_ID) %>%
  mutate(Genome_ARG_Count = sum(Present, na.rm = TRUE)) %>%
  ungroup()




gene_class_lookup <- ARG_entero_plasmid %>%
  select(gene_id, Gene_class) %>%
  distinct()  # in case there are duplicates

Plasmid_Source <- Plasmid_Source %>%
  left_join(gene_class_lookup, by = "gene_id")

Plasmid_Source <- Plasmid_Source %>%
  group_by(Genome_ID) %>%                                  # group by genome
  mutate(BETA_LACTAM_present = if_else(
    any(Gene_class == "BETA-LACTAM"),                     # check if any gene_class is BETA-LACTAM
    "Y", 
    "N"
  )) %>%
  ungroup()


Plasmid_Source_Beta <- Plasmid_Source %>%
  filter(BETA_LACTAM_present == "Y")

Plasmid_Source_Beta <- Plasmid_Source_Beta %>%
  mutate(Gene_class = factor(Gene_class, 
                             levels = c("BETA-LACTAM", 
                                        setdiff(unique(Gene_class), "BETA-LACTAM"))))

Plasmid_Source_Beta <- Plasmid_Source_Beta %>%
  mutate(gene_id = case_when(
    gene_id %in% c("oqxA", "oqxB") ~ "oqxAB",
    gene_id %in% c("dfrA14", "dfrA17") ~ "dfrA",
    TRUE ~ gene_id  # keep all other genes unchanged
  ))
unique(Plasmid_Source_Beta$gene_id)

Plasmid_Source_Beta <- Plasmid_Source_Beta %>%
  mutate(gene_id = factor(gene_id, levels = c(
    # BETA-LACTAM first, alphabetically
    "blaCTX-M-15", "blaCTX-M-27", "blaEC", "blaLEN-16", "blaOXA-1", "blaSHV-27", "blaTEM-1",
    # Next class (example: aminoglycosides), alphabetically
    "aac(6')-Ib-cr5", "aadA5", "aph(3'')-Ib", "aph(6)-Id",
    # Sulfonamides
    "sul1", "sul2",
    # Chloramphenicol
    "catA1",
    # Tetracyclines
    "tet(A)", "tet(B)",
    # Efflux / others
    "kdeA", "mdtM", "oqxAB",
    # Additional classes
    "dfrA", "mph(A)", "qnrB1"
  )))


Plasmid_Source_Beta <- Plasmid_Source_Beta %>%
  mutate(`Source2` = case_when(
    grepl("pContig", Source_Genome) ~ "Contig",
    grepl("pMAG", Source_Genome) ~ "MAG",
    TRUE ~ "Other"  # optional, for anything else
  ))

Plasmid_Source_Beta <- Plasmid_Source_Beta %>%
  arrange(Genus, Genome_ID) %>%
  mutate(Genome_ID = factor(Genome_ID, levels = unique(Genome_ID)))

unique(Plasmid_Source_Beta$Genome_ID)

Genome_labels <- c(
  "CD01_bin_7" = "CD01 (bin_7)",
  "CD05_bin_7" = "CD05 (bin_7)",
  "CD07_bin_9"     = "CD07 (bin_9)",
  "CD14_bin_10"    = "CD14 (bin_10)",
  "CD15_bin_6"     = "CD15 (bin_6)",
  "CD23_bin_3"     = "CD23 (bin_3)",
  "CD24_bin_8"     = "CD24 (bin_8)",
  "CD26_bin_4"     = "CD26 (bin_4)", #This originally had bin_24 which in label name which is incorrect
  "CD31_bin_6"     = "CD31 (bin_6)",
  "CD33_bin_5"     = "CD33 (bin_5)",
  "CD04_bin_8" = "CD04 (bin_8)",
  "CD29_bin_10"    = "CD29 (bin_10)",
  "CD34_bin_7"     = "CD34 (bin_7)",
  "CD02_bin_8" = "CD02 (bin_8)"
)


Fig2A_1 <-ggplot(Plasmid_Source_Beta, aes(x = gene_id, y = Genome_ID)) +
  geom_tile(fill = "grey70", color = "black", linewidth = 0.2) +  # simple grey fill
  theme_classic() +
  scale_y_discrete(labels = Genome_labels) +  
  labs(x = "Gene Group", y = "Genome") +  # no legend needed
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  ) + coord_fixed(ratio = 1)



#####Code to plot Enterobacteriaceae VFs in MAGs with plasmid borne beta-lactamases
vf_entero_MDR <- vf_entero %>%
  filter(Genome_ID %in% Plasmid_Source_Beta$Genome_ID)

vf_entero_MDR <- vf_entero_MDR %>%
  mutate(Genome_ID = factor(Genome_ID,
                            levels = levels(Plasmid_Source_Beta$Genome_ID)))

# Get all Genome_ID levels
all_genomes <- levels(Plasmid_Source_Beta$Genome_ID)

# Identify missing Genome_IDs
missing_genomes <- setdiff(all_genomes, unique(vf_entero_MDR$Genome_ID))

# Create one blank row per missing genome
blank_rows <- tibble(
  Genome_ID = factor(missing_genomes, levels = all_genomes),
  gene_id   = NA
)

# Add blank rows to your data
vf_entero_MDR_padded <- bind_rows(
  vf_entero_MDR %>%
    mutate(Genome_ID = factor(Genome_ID, levels = all_genomes)),
  blank_rows
)

unique(vf_entero_MDR_padded$gene_id)


vf_entero_MDR_padded <- vf_entero_MDR_padded %>%
  mutate(
    VF_group = case_when(
      gene_id %in% c("cnf1", "hlyA-alpha", "sat", "vactox", "senB") ~ "Toxins",
      gene_id %in% c("iroB","iroC","iroD","iroE","iroN") ~ "Salmochelin",
      gene_id %in% c("iucA","iucB","iucC","iucD","iutA") ~ "Aerobactin",
      gene_id %in% c("ybtP","ybtQ") ~ "Yersiniabactin",
      gene_id %in% c("papC","papF","papG-II","papG-III","papH", "papE") ~ "P fimbriae",
      gene_id %in% c("sfaF","sfaS") ~ "S fimbriae",
      gene_id %in% c("lpfA") ~ "long polar fimbriae",
      gene_id %in% c("eilA","espX1","iha","fdeC", "focG") ~ "Other adhesins",
      gene_id %in% c("iss","sslE","ibeA","mchF") ~ "Other",
      gene_id %in% c("cvaC") ~ "Colicin"
    )
  )

vf_entero_MDR_padded <- vf_entero_MDR_padded %>%
  mutate(VF_group = factor(VF_group, levels = c(
    "Toxins",
    "P fimbriae",
    "S fimbriae",
    "long polar fimbriae",
    "Other adhesins",
    "Salmochelin",
    "Aerobactin",
    "Yersiniabactin",
    "Colicin",
    "Other"
  )))


vf_plot_data <- vf_entero_MDR_padded %>%
  group_by(Genome_ID, VF_group) %>%
  summarise(
    total_genes = n(),
    plasmid_genes = sum(Contig_type == "plasmid", na.rm = TRUE),
    plasmid_fraction = plasmid_genes / total_genes
  ) %>%
  ungroup()


#Note I needed to add the "NA" x axis to alin fgure with other heatmap. This is actually a placeholder for the strains without VFs but have plasmid beta-lactamases. This row is removed in graphics editor 
Fig2A_2 <- ggplot(vf_plot_data, aes(x = VF_group, y = Genome_ID)) +
  geom_tile(aes(fill = plasmid_fraction), color = "black", linewidth = 0.2) +
  geom_text(aes(label = ifelse(total_genes > 1, total_genes, "")), size = 3) +
  scale_fill_gradient(low = "white", high = "grey70", name = "% plasmid") +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  ) +
  labs(x = "VF Group", y = NULL) + 
  theme(legend.position = "top")

#####Combine and Plot
Fig2A <- Fig2A_1 + Fig2A_2
Fig2A
