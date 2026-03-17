###Author: Andrew J Sommer
###Date Modified: 3/17/2026
###Description: Code for Analysis of ARG Accumulation by Healthcare Burden 
library(dplyr)
library(tidyverse)
library(ggplot2)
library(rstatix)

setwd("~/Work/UMN-Postdoc/Projects/2025/2025-PhaseGenomics/Github-Upload/")

#####Load in metadata
sample_meta <- read.csv("Import/Metadata/Sample-Metadata.csv")
head(sample_meta)


####ARG Summary 
arg_cleaned <- read.csv("Import/Data/AMR_VF/arg-table-core.csv")
ARG_entero <- arg_cleaned %>%  filter(Family == "Enterobacteriaceae")
summary_entero_arg <- ARG_entero %>%
  complete(ID = sample_meta$ID) %>%  
  group_by(ID) %>%
  summarise(Sum_Contig_TPM = sum(Contig_TPM, na.rm = TRUE)) %>%
  mutate(log1p_Sum = log1p(Sum_Contig_TPM)) %>%
  left_join(sample_meta, by = "ID")
summary_entero_arg <- summary_entero_arg  %>%
  rename(ARG_TPM = Sum_Contig_TPM) %>%
  rename(ARG_log1p = log1p_Sum) %>%
  select(-treatment, -alias)

### VF Summary
vf_table_cleaned <- read.csv(file = "Import/Data/AMR_VF/virulence-gene-table.csv") #cov.adj > 90%
vf_entero <- vf_table_cleaned %>%  filter(Family == "Enterobacteriaceae")
summary_entero_vf <- vf_entero %>%
  complete(ID = sample_meta$ID) %>%  
  group_by(ID) %>%
  summarise(Sum_Contig_TPM = sum(Contig_TPM, na.rm = TRUE)) %>%
  mutate(log1p_Sum = log1p(Sum_Contig_TPM)) %>%
  left_join(sample_meta, by = "ID")
summary_entero_vf <- summary_entero_vf  %>%
  rename(VF_TPM = Sum_Contig_TPM) %>%
  rename(VF_log1p = log1p_Sum) 


### Combine into a scatterplot
scatter_df <- summary_entero_arg %>%
  left_join(summary_entero_vf, by = c("ID"))




###############################################Extended Fig 6: rCDI Patient Burden Analysis#############################################

#### Combine VF/ARG summary with rCDI Metadata
rcdi_patient_meta <- read.csv("Import/Metadata/Patient-Metadata-rCDI-Key.csv")
rcdi_burden <- rcdi_patient_meta %>% left_join(scatter_df, by = c("ID"))
rcdi_burden <- rcdi_burden %>% mutate(burden = factor(burden, levels = c("Low", "Medium",  "High")))


#### Boxplots of VF and ARG Abundance 
kruskal.test(VF_log1p ~ burden, data = rcdi_burden)
kruskal.test(ARG_log1p ~ burden, data = rcdi_burden)

CDI_VF_Boxplot <- ggplot(rcdi_burden, aes(x = burden, y = VF_log1p, fill = burden)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  labs(
    x = "Patient Burden",
    y = "VF abundance\nlog1p(TPM)") +
  theme_classic() +
  scale_fill_manual(
    values = c(
      "Low"    = "#1b9e77",  # green
      "Medium" = "#d95f02",  # orange
      "High"   = "red")) +
  theme(legend.position = "none")

CDI_ARG_Boxplot <- ggplot(rcdi_burden, aes(x = burden, y = ARG_log1p, fill = burden)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  labs(
    x = "Burden category",
    y = "ARG abundance\nlog1p(TPM)") +
  theme_classic() +
  scale_fill_manual(
    values = c(
      "Low"    = "#1b9e77",  # green
      "Medium" = "#d95f02",  # orange
      "High"   = "red")) +
  theme(legend.position = "none")



###rCDI Saccterplots. Note that this plots only patients with the available metadata
rcdi_burden$Duration_CDI <- ifelse(rcdi_burden$Duration_CDI == "> 3 years", "36",
                                   ifelse(rcdi_burden$Duration_CDI == "> 12 months", "12",
                                          rcdi_burden$Duration_CDI))
rcdi_burden$Duration_CDI_numeric <- as.numeric(rcdi_burden$Duration_CDI)

rcdi_burden$Age <- as.numeric(rcdi_burden$Age)
cor.test(rcdi_burden$Duration_CDI_numeric, rcdi_burden$VF_log1p, method = "spearman")

CDI_Scatter_A <- ggplot(rcdi_burden, aes(x = Duration_CDI_numeric, y =VF_log1p,fill = burden)) +
  geom_jitter(
    shape = 21,       # circle with border
    width = 0.1,
    height = 0.1,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    x = "Duration rCDI Treatment\n(months)",
    y = "Enterobacteriaceae\nVFs log1p(TPM)",
    fill = "Burden",
    subtitle = "p = 0.88, rho=-0.28"
  ) +
  theme_classic() +
  scale_fill_manual(
    values = c(
      "Low"    = "#1b9e77",  # green
      "Medium" = "#d95f02",  # orange
      "High"   = "red")) +
  theme(legend.position = "none")



cor.test(rcdi_burden$Duration_CDI_numeric, rcdi_burden$ARG_log1p, method = "spearman")
CDI_Scatter_B <- ggplot(rcdi_burden, aes(x = Duration_CDI_numeric, y =ARG_log1p,fill = burden)) +
  geom_jitter(
    shape = 21,       # circle with border
    width = 0.1,
    height = 0.1,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    x = "Duration rCDI Treatment\n(months)",
    y = "Enterobacteriaceae\nARGs log1p(TPM)",
    fill = "Burden",
    subtitle = "p = 0.39, rho=0.153"
  ) +
  theme_classic() +
  scale_fill_manual(
    values = c(
      "Low"    = "#1b9e77",  # green
      "Medium" = "#d95f02",  # orange
      "High"   = "red"))+
  theme(legend.position = "none")




### Age of Patient
cor.test(rcdi_burden$Age, rcdi_burden$VF_log1p, method = "spearman")
CDI_Scatter_C <- ggplot(rcdi_burden, aes(x = Age, y =VF_log1p,fill = burden)) +
  geom_jitter(
    shape = 21,       # circle with border
    width = 0.1,
    height = 0.1,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    x = "Age",
    y = "Enterobacteriaceae\nVFs log1p(TPM)",
    fill = "Burden",
    subtitle = "p = 0.14, rho=0.26"
  ) +
  theme_classic() +
  scale_fill_manual(
    values = c(
      "Low"    = "#1b9e77",  # green
      "Medium" = "#d95f02",  # orange
      "High"   = "red")) +
  theme(legend.position = "none")


cor.test(rcdi_burden$Age, rcdi_burden$ARG_log1p, method = "spearman")
CDI_Scatter_D <- ggplot(rcdi_burden, aes(x = Age, y =ARG_log1p,fill = burden)) +
  geom_jitter(
    shape = 21,       # circle with border
    width = 0.1,
    height = 0.1,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    x = "Age",
    y = "Enterobacteriaceae\nARGs log1p(TPM)",
    fill = "Burden",
    subtitle = "p = 0.08, rho=0.30"
  ) +
  theme_classic() +
  scale_fill_manual(
    values = c(
      "Low"    = "#1b9e77",  # green
      "Medium" = "#d95f02",  # orange
      "High"   = "red")) +
  theme(legend.position = "none")

library(patchwork)
combined_plot <- (CDI_VF_Boxplot | CDI_ARG_Boxplot) / (CDI_Scatter_A | CDI_Scatter_B) / (CDI_Scatter_C | CDI_Scatter_D) + 
  plot_annotation(
    tag_levels = 'A'                # automatically label plots as A, B, C, D
  )

combined_plot


###############################################Extended Fig 7: Liver Patient Burden Analysis#############################################


###Load in Liver Patient Metdata
liver_key <- read.csv("Import/Metadata/Patient-Metadata-Liver-Key.csv")
liver_key <- liver_key %>% dplyr::select(ID,Age,Rifaximin)
liver_meta <- left_join(liver_key,sample_meta, by ="ID")

liver_burden <- liver_meta %>% left_join(scatter_df, by = c("ID", "alias", "treatment"))
colnames(liver_burden)

color_code <- c(
  "Compensated"    = "#f0e548ff",
  "Decompensated"  = "#e69696ff")


CDI_ARG_Boxplot <- ggplot(liver_burden, 
                  aes(x = Rifaximin, y = ARG_log1p, fill = treatment)) +
  geom_boxplot(aes(fill = Rifaximin), alpha = 0.1, outlier.shape = NA) +  
  geom_jitter(width = 0.2, size = 4, shape = 21, color = "black") +
  labs(y = "Enterobacteriaceae ARGs\nlog1p(TPM)", x = "",
       subtitle = "p = 0.240") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 11),
    legend.position = "none"
  ) + 
  scale_x_discrete(labels = c(
    "N" = "No Rifaximin Treatment",
    "Y" = "Rifaximin Treatment"))  +
  scale_fill_manual(values = color_code)
CDI_ARG_Boxplot

CDI_VF_Boxplot <- ggplot(liver_burden, 
                 aes(x = Rifaximin, y = VF_log1p, fill = treatment)) +
  geom_boxplot(aes(fill = Rifaximin), alpha = 0.1, outlier.shape = NA) +  
  geom_jitter(width = 0.2, size = 4, shape = 21, color = "black") +
  labs(y = "Enterobacteriaceae VFs\nlog1p(TPM)", x = "",
       subtitle = "p = 0.328") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 11),
    legend.position = "none"
  ) + 
  scale_x_discrete(labels = c(
    "N" = "No Rifaximin Treatment",
    "Y" = "Rifaximin Treatment"))  +
  scale_fill_manual(values = color_code)
CDI_VF_Boxplot



wilcox.test(VF_log1p ~ treatment, data = liver_burden) # p =0.24
wilcox.test(ARG_log1p ~ treatment, data = liver_burden) # p =0.24



CDI_VF_Scatter <- ggplot(liver_burden, aes(x = Age, y = VF_log1p, fill = treatment)) +
  geom_jitter(
    shape = 21,       # circle with border
    color = "black",  # border color
    width = 0.1,
    height = 0.1,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    x = "Patient Age (Years)",
    y = "Enterobacteriaceae\nVFs log1p(TPM)",
    fill = "Burden",
    subtitle = "p = 0.27, rho=0.-267"
  ) +
  theme_classic() +
  scale_fill_manual(values = color_code) +
  theme(legend.position = "none")


CDI_ARG_Scatter <- ggplot(liver_burden, aes(x = Age, y = ARG_log1p, fill = treatment)) +
  geom_jitter(
    shape = 21,       # circle with border
    color = "black",  # border color
    width = 0.1,
    height = 0.1,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    x = "Patient Age (Years)",
    y = "Enterobacteriaceae\nARGs log1p(TPM)",
    fill = "Burden",
    subtitle = "p = 0.996, rho=0.001"
  ) +
  theme_classic() +
  scale_fill_manual(values = color_code) +
  theme(legend.position = "none")

cor.test(liver_burden$Age, liver_burden$VF_log1p, method = "spearman")
cor.test(liver_burden$Age, liver_burden$ARG_log1p, method = "spearman")


combined_plot <- (CDI_VF_Boxplot | CDI_ARG_Boxplot) / (CDI_VF_Scatter | CDI_ARG_Scatter)+ 
  plot_annotation(
    tag_levels = 'A'                # automatically label plots as A, B, C, D
  )

combined_plot

###############################################Extended Fig 8: Combined Rifaximin Analysis#############################################

###Set Colors
color_code <- c(
  "Decompensated"  = "#e69696ff",
  "rCDI"           = "#2b6e88ff")

###ARG/VF
setwd("~/Work/UMN-Postdoc/Projects/2025/2025-PhaseGenomics/Analysis/Final/Import-New/AMR_VF")
ARG_VF_Summary <- read.csv("ARG_VF_Summary.csv")



ARG_VF_Summary_meta <- ARG_VF_Summary %>% left_join(sample_meta, by ="ID")
ARG_VF_Summary_meta_select <- ARG_VF_Summary_meta %>% filter(treatment != "Compensated")
ARG_VF_Summary_meta_select <- ARG_VF_Summary_meta_select %>% filter(treatment != "Healthy")


###################Enterobacteriaceae VFs
wilcox.test(Enterobacteriacea_VF_lop1p ~ Rifamixin, data = ARG_VF_Summary_meta_select)

plotVF <- ggplot(ARG_VF_Summary_meta_select, 
                 aes(x = Rifamixin, y = Enterobacteriacea_VF_lop1p, fill = treatment)) +
  geom_boxplot(aes(fill = Rifamixin), alpha = 0.1, outlier.shape = NA) +  
  geom_jitter(width = 0.2, size = 4, shape = 21, color = "black") +
  labs(y = "Enterobacteriaceae VFs\nlog1p(TPM)", x = "",
       subtitle = "p = 0.274") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 11),
    legend.position = "none"
  ) + 
  scale_x_discrete(labels = c(
    "N" = "No Rifaximin Treatment",
    "Y" = "Prior Rifaximin Treatment"))  +
  scale_fill_manual(values = color_code)
plotVF


###################Enterobacteriaceae ARGs
wilcox.test(Enterobacteriaceae_ARG_lop1p ~ Rifamixin, data = ARG_VF_Summary_meta_select)
plotARG <- ggplot(ARG_VF_Summary_meta_select, 
                  aes(x = Rifamixin, y = Enterobacteriaceae_ARG_lop1p, fill = treatment)) +
  geom_boxplot(aes(fill = Rifamixin), alpha = 0.1, outlier.shape = NA) +  
  geom_jitter(width = 0.2, size = 4, shape = 21, color = "black") +
  labs(y = "Enterobacteriaceae ARGs\nlog1p(TPM)", x = "",
       subtitle = "p = 0.691") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 11),
    legend.position = "none"
  ) + 
  scale_x_discrete(labels = c(
    "N" = "No Rifaximin Treatment",
    "Y" = "Prior Rifaximin Treatment"))  +
  scale_fill_manual(values = color_code)
plotARG



############Vancomycin Resistance
wilcox.test(Vancomycin_ARGs_log1p ~ Rifamixin, data = ARG_VF_Summary_meta_select)
plotVanco <- ggplot(ARG_VF_Summary_meta_select, 
                    aes(x = Rifamixin, y = Vancomycin_ARGs_log1p, fill = treatment)) +
  geom_boxplot(aes(fill = Rifamixin), alpha = 0.1, outlier.shape = NA) +  
  geom_jitter(width = 0.2, size = 4, shape = 21, color = "black") +
  labs(y = "Vancomycin ARGs\nlog1p(TPM)", x = "",
       subtitle = "p = 0.651") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 11),
    legend.position = "none"
  ) + 
  scale_x_discrete(labels = c(
    "N" = "No Rifaximin Treatment",
    "Y" = "Prior Rifaximin Treatment"))  +
  scale_fill_manual(values = color_code)
plotVanco


library(patchwork)
final_plot <-
  (plotVF +plotARG  + plotVanco) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
final_plot












