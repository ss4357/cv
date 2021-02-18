library(tidyverse)
library(readxl)
library(ggfortify)
library(ggplot2)
library(EnhancedVolcano)

project_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

#ms2 validations---------------------------------------------------------------------------------
valid_lp <- read_xlsx(paste0(project_dir, "/data/2020-05-12/lipidomics_pos_ms2_validation_results.xlsx")) %>%
  filter(!is.na(name))
valid_ln <- read_xlsx(paste0(project_dir, "/data/2020-05-12/lipidomics_neg_ms2_validation_results.xlsx")) %>%
  filter(!is.na(name))
valid_mp <- read_xlsx(paste0(project_dir, "/data/2020-05-12/metabolomics_pos_ms2_validation_results.xlsx")) %>%
  filter(!is.na(name))
valid_mn <- read_xlsx(paste0(project_dir, "/data/2020-05-12/metabolomics_neg_ms2_validation_results.xlsx")) %>%
  filter(!is.na(name))

# PCA --------------------------------------------------------------------------------

pcaplots <- function(mode, comparison){
  comp_table <- read_csv(paste0(project_dir, "/data/", mode, "/", comparison), col_names = FALSE, skip = 2)
  nums <- comp_table[,-1]
  pca_temp <- prcomp(nums, scale. = TRUE)
  plot_temp <- autoplot(pca_temp, colour = "coral2")
  return(plot_temp)
}

rad_v_sham <- "1. radiation_vs_sham.csv"
veh_v_sham <- "2. vehicle_sd4_vs_sham.csv"
exrad1_v_veh <- "3. exrad1_sd4_vs_vehicle_sd4.csv"
exrad2_v_veh <- "4. exrad2_sd4_vs_vehicle_sd4.csv"

ln <- "exrad_serum_lipidomics_NEG"
lp <- "exrad_serum_lipidomics_POS"
mn <- "exrad_serum_metabolomics_NEG"
mp <- "exrad_serum_metabolomics_POS"

jpeg(paste0(project_dir, "/ln_rad_v_sham.jpg"))
pcaplots(ln, rad_v_sham)
dev.off()
pcaplots(ln, veh_v_sham)
pcaplots(ln, exrad1_v_veh)
pcaplots(ln, exrad2_v_veh)

pcaplots(lp, rad_v_sham)
pcaplots(lp, veh_v_sham)
pcaplots(lp, exrad1_v_veh)
pcaplots(lp, exrad2_v_veh)

pcaplots(mn, rad_v_sham)
pcaplots(mn, veh_v_sham)
pcaplots(mn, exrad1_v_veh)
pcaplots(mn, exrad2_v_veh)

pcaplots(mp, rad_v_sham)
pcaplots(mp, veh_v_sham)
pcaplots(mp, exrad1_v_veh)
pcaplots(mp, exrad2_v_veh)

# Volcano -----------------------------------------------------------------------------------------
# The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
stats_lp <- read_xlsx(paste0(project_dir, "/data/exrad_serum_lipidomics_POS/stats_all.xlsx"), 
                      col_names = c("index", "name", 
                                    "p_value_x", "FDR_x", "foldchange_x", "log2FC_x",
                                    "p_value_y", "FDR_y", "foldchange_y", "log2FC_y", 
                                    "p_value_x_x", "FDR_x_x", "foldchange_x_x", "log2FC_x_x",
                                    "p_value_y_y", "FDR_y_y", "foldchange_y_y", "log2FC_y_y"), skip = 1)
stats_ln <- read_xlsx(paste0(project_dir, "/data/exrad_serum_lipidomics_NEG/stats_all.xlsx"), 
                      col_names = c("index", "name", 
                                    "p_value_x", "FDR_x", "foldchange_x", "log2FC_x",
                                    "p_value_y", "FDR_y", "foldchange_y", "log2FC_y", 
                                    "p_value_x_x", "FDR_x_x", "foldchange_x_x", "log2FC_x_x",
                                    "p_value_y_y", "FDR_y_y", "foldchange_y_y", "log2FC_y_y"), skip = 1)
stats_mp <- read_xlsx(paste0(project_dir, "/data/exrad_serum_metabolomics_POS/stats_all.xlsx"), 
                      col_names = c("index", "name", 
                                    "p_value_x", "FDR_x", "foldchange_x", "log2FC_x",
                                    "p_value_y", "FDR_y", "foldchange_y", "log2FC_y", 
                                    "p_value_x_x", "FDR_x_x", "foldchange_x_x", "log2FC_x_x",
                                    "p_value_y_y", "FDR_y_y", "foldchange_y_y", "log2FC_y_y"), skip = 1)
stats_mn <- read_xlsx(paste0(project_dir, "/data/exrad_serum_metabolomics_NEG/stats_all.xlsx"), 
                      col_names = c("index", "name", 
                                    "p_value_x", "FDR_x", "foldchange_x", "log2FC_x",
                                    "p_value_y", "FDR_y", "foldchange_y", "log2FC_y", 
                                    "p_value_x_x", "FDR_x_x", "foldchange_x_x", "log2FC_x_x",
                                    "p_value_y_y", "FDR_y_y", "foldchange_y_y", "log2FC_y_y"), skip = 1)

valid_stats_lp <- stats_lp %>% left_join(valid_lp %>% dplyr::select(3, 4) 
                                         %>% dplyr::rename(chem_name=name), by = c("name" = "precursor_mz_rt"))
valid_stats_ln <- stats_ln %>% left_join(valid_ln %>% dplyr::select(3, 4) 
                                         %>% dplyr::rename(chem_name=name), by = c("name" = "precursor_mz_rt"))
valid_stats_mp <- stats_mp %>% left_join(valid_mp %>% dplyr::select(3, 4) 
                                         %>% dplyr::rename(chem_name=name), by = c("name" = "precursor_mz_rt"))
valid_stats_mn <- stats_mn %>% left_join(valid_mn %>% dplyr::select(3, 4) 
                                         %>% dplyr::rename(chem_name=name), by = c("name" = "precursor_mz_rt"))


all_stats <- bind_rows(valid_stats_lp, valid_stats_ln, valid_stats_mp, valid_stats_mn)

# Rad v sham volcano
all_stats %>% dplyr::select(2,3,6)  %>%
  EnhancedVolcano(lab = all_stats$chem_name, 
                  x = 'log2FC_x', 
                  y = "p_value_x",
                  title = "Radiation Verus Sham")

# Vehicle SD4 vs Sham volcano
all_stats %>% dplyr::select(2,7,10)  %>%
  EnhancedVolcano(lab = all_stats$name, 
                  x = 'log2(FC).y', 
                  y = "p.value.y",
                  title = "Vehicle SD4 Versus Sham")

#EXRADI SD4 Versus Vehicle SD4 volcano
all_stats %>% dplyr::select(2, 11, 14)  %>%
  EnhancedVolcano(lab = all_stats$chem_name, 
                  x = 'log2FC_x_x', 
                  y = "p_value_x_x")

#EXRADI SD4 Versus Vehicle SD4 volcano
all_stats %>% dplyr::select(2, 15, 18)  %>%
  EnhancedVolcano(lab = all_stats$chem_name, 
                  x = 'log2(FC).y.y', 
                  y = "p.value.y.y",
                  title = "EXRADII SD4 Versus Vehicle SD4")


# Biomarker identification ------------------------------------------------------------------------------
valid <- bind_rows(valid_ln, valid_lp, valid_mn, valid_mp)

#match validations to stats, by mz_rt
valid <- valid %>%
  dplyr::select(mode, precursor_mz_rt, name) %>%
  left_join(all_stats, by = c("precursor_mz_rt" = "name")) %>%
  dplyr::select(1:3, 8, 12, 16, 20)

#drop column with NA
valid <- valid[-21,]

pheatmap(as.matrix(valid[,4:7]), labels_col = c("rad_v_sham", "vehicle_sd4_v_sham", "exrad1_v_vehicle", "exrad2_v_vehicle"))

pos_rad_v_sham <- filter(valid, valid[,4] > 0)
neg_rad_v_sham <- filter(valid, valid[,4] < 0)
pos_reg_rad_v_sham <- pos_rad_v_sham[,3]
neg_reg_rad_v_sham <- neg_rad_v_sham[,3]

write_tsv(pos_reg_rad_v_sham, paste0(project_dir, 
                                     "/pos_reg_rad_v_sham.tsv"))
write_tsv(neg_reg_rad_v_sham, paste0(project_dir, 
                                     "/neg_reg_rad_v_sham.tsv"))
