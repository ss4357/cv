library(tidyverse)
library(readxl)

project_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

#function: filter function for p value < 0.05 -----------------------------------------------
filtering <- function(all){
  rad_filtered <- filter(all, p_value < 0.05)
  rad_filtered <- arrange(rad_filtered, pathway)
  
  rad_filt_combined_str <- rad_filtered %>% 
    group_by(pathway) %>%
    mutate(p_value = paste0(p_value),
           empirical_compounds = paste0(empirical_compounds, collapse = ""), 
           feature_id = paste0(feature_id, collapse = ""),
           feature_name = paste0(feature_name, collapse = ""))
  
  filt_sum_rad <- rad_filt_combined_str %>%
    group_by(pathway, empirical_compounds, feature_id, feature_name) %>%  
    summarise(overlap = sum(overlap), 
              pathway_size = sum(pathway_size),
              p_value = min(p_value))
  return(filt_sum_rad)
}

#function: get empriical ID function ---------------------------------------------------------------
get_eid <- function(mode_data, emp_data){
  mode_data <- separate_rows(mode_data, empirical_compounds, sep = c(","))
  lp_filt <- filter(mode_data, p_value < 0.05, overlap > 0, pathway_size > 0)
  lp_filt <- left_join(select(lp_filt, pathway, empirical_compounds), emp_data, by = c("empirical_compounds" = "EID"))
  lp_filt <- separate_rows(lp_filt, massfeature_rows)
  return(lp_filt)
}

#function: match pathway comppounds with MS2 validations ----------------------------------------------
match <- function(link, mode){
  #get pathway analysis tables
  read_mode <- read_tsv(paste0(project_dir, link, "mcg_pathwayanalysis_pathway.tsv"),
                        col_names = c("pathway", "overlap", "pathway_size", "p_value", 
                                      "empirical_compounds", "feature_id", "feature_name"), skip = 1)
  
  #get emprical compounds
  emp <- read_tsv(paste0(project_dir, link, "ListOfEmpiricalCompounds.tsv"))
  
  #read user input data, which has mz-rt
  input_data <- read_tsv(paste0(project_dir, link, "userInputData.txt"))

  #add empirical compounds EID/row num to pathway table
  filt <- get_eid(read_mode, emp) %>%
    left_join(select(input_data, massfeature_rows, CompoundID_from_user), by = c("massfeature_rows")) %>%
    left_join(select(get(paste0("valid_", mode)), mode, precursor_mz_rt, name, synonyms, adduct, level),
              by = c("CompoundID_from_user" = "precursor_mz_rt")) %>%
    distinct()
  return(filt)
}

#function: read & combine tables to get all pathway compounds --------------------------------------------------------------------------
combine_tables <- function(lp, ln, mp, mn){
  lipid_pos <- read_tsv(paste0(project_dir, lp, "mcg_pathwayanalysis_pathway.tsv"),
                        col_names = c("pathway", "overlap", "pathway_size", "p_value", "empirical_compounds", "feature_id", "feature_name"), skip = 1)
  lipid_neg <- read_tsv(paste0(project_dir, ln, "mcg_pathwayanalysis_pathway.tsv"),
                        col_names = c("pathway", "overlap", "pathway_size", "p_value", "empirical_compounds", "feature_id", "feature_name"), skip = 1)
  met_pos <- read_tsv(paste0(project_dir, mp, "mcg_pathwayanalysis_pathway.tsv"),
                      col_names = c("pathway", "overlap", "pathway_size", "p_value", "empirical_compounds", "feature_id", "feature_name"), skip = 1)
  met_neg <- read_tsv(paste0(project_dir, mn, "mcg_pathwayanalysis_pathway.tsv"),
                      col_names = c("pathway", "overlap", "pathway_size", "p_value", "empirical_compounds", "feature_id", "feature_name"), skip = 1)
  temp_table <-  bind_rows(lipid_pos, lipid_neg, met_neg, met_pos)
  return(temp_table)
}


#function: add not validated column if there is no MS2 value -------------------------------------
match_ms2_pathway <- function(table){
  table <- table %>%  
    mutate(match = ifelse(is.na(name), "Not Validated", "")) %>%
    select(CompoundID_from_user, pathway_compound_names = compound_names, ms2_compound_names = name,match, pathway_compounds = compounds, pathway, ms2_synonyms = synonyms, )
  return(table)  
}

#ms2 validations---------------------------------------------------------------------------------
valid_lp <- read_xlsx(paste0(project_dir, "/data/2020-05-12/lipidomics_pos_ms2_validation_results.xlsx")) %>%
  filter(!is.na(name))
valid_ln <- read_xlsx(paste0(project_dir, "/data/2020-05-12/lipidomics_neg_ms2_validation_results.xlsx")) %>%
  filter(!is.na(name))
valid_mp <- read_xlsx(paste0(project_dir, "/data/2020-05-12/metabolomics_pos_ms2_validation_results.xlsx")) %>%
  filter(!is.na(name))
valid_mn <- read_xlsx(paste0(project_dir, "/data/2020-05-12/metabolomics_neg_ms2_validation_results.xlsx")) %>%
  filter(!is.na(name))

# Radiation vs Sham --------------------------------------------------------------

lp <- "/data/exrad_serum_lipidomics_POS/1. radiation_vs_sham/1583349765.12.pathway/tables/"
ln <- "/data/exrad_serum_lipidomics_NEG/1. radiation_vs_sham/1583349669.51.pathway/tables/"
mp <- "/data/exrad_serum_metabolomics_POS/1. radiation_vs_sham/1583360360.97.pathway/tables/"
mn <- "/data/exrad_serum_metabolomics_NEG/1. radiation_vs_sham/1583360271.07.pathway/tables/"

match_lp <- match(lp, x<-"lp")
match_ln <- match(ln, x<-"ln")
match_mp <- match(mp, x<-"mp")
match_mn <- match(mn, x<-"mn")

compound_comp_rad_v_sham <- bind_rows(match_lp, match_ln, match_mp, match_mn)

#rad v sham compounds filtered

all <- combine_tables(lp, ln, mp, mn)
filtered_rad_v_sham <- filtering(all)

# Vehicle vs Rad --------------------------------------------------------------
lp <- "/data/exrad_serum_lipidomics_POS/2. vehicle_sd4_vs_sham/1583349788.37.pathway/tables/"
ln <- "/data/exrad_serum_lipidomics_NEG/2. vehicle_sd4_vs_sham/1583349697.19.pathway/tables/"
mp <- "/data/exrad_serum_metabolomics_POS/2. vehicle_sd4_vs_sham/1583360392.09.pathway/tables/"
mn <- "/data/exrad_serum_metabolomics_NEG/2. vehicle_sd4_vs_sham/1583360295.83.pathway/tables/"

match_lp <- match(lp, x<-"lp")
match_ln <- match(ln, x<-"ln")
match_mp <- match(mp, x<-"mp")
match_mn <- match(mn, x<-"mn")

compound_comp_vehicle_v_rad <- bind_rows(match_lp, match_ln, match_mp, match_mn)

all <- combine_tables(lp, ln, mp, mn)
filtered_veh_v_rad <- filtering(all)


# ExRADI vs Rad --------------------------------------------------------------
lp <- "/data/exrad_serum_lipidomics_POS/3. exrad1_sd4_vs_vehicle_sd4/1583349811.16.pathway/tables/"
ln <- "/data/exrad_serum_lipidomics_NEG/3. exrad1_sd4_vs_vehicle_sd4/1583349724.11.pathway/tables/"
mp <- "/data/exrad_serum_metabolomics_POS/3. exrad1_sd4_vs_vehicle_sd4/1583360421.95.pathway/tables/"
mn <- "/data/exrad_serum_metabolomics_NEG/3. exrad1_sd4_vs_vehicle_sd4/1583360319.45.pathway/tables/"

match_lp <- match(lp, x<-"lp")
match_ln <- match(ln, x<-"ln")
match_mp <- match(mp, x<-"mp")
match_mn <- match(mn, x<-"mn")

compound_comp_exradi <- bind_rows(match_lp, match_ln, match_mp, match_mn)

all <- combine_tables(lp, ln, mp, mn)
filtered_exradi_v_rad <- filtering(all)

# ExRADII vs Rad --------------------------------------------------------------
lp <- "/data/exrad_serum_lipidomics_POS/4. exrad2_sd4_vs_vehicle_sd4/1583349824.94.pathway/tables/"
ln <- "/data/exrad_serum_lipidomics_NEG/4. exrad2_sd4_vs_vehicle_sd4/1583349745.33.pathway/tables/"
mp <- "/data/exrad_serum_metabolomics_POS/4. exrad2_sd4_vs_vehicle_sd4/1583360446.38.pathway/tables/"
mn <- "/data/exrad_serum_metabolomics_NEG/4. exrad2_sd4_vs_vehicle_sd4/1583360340.23.pathway/tables/"

match_lp <- match(lp, x<-"lp")
match_ln <- match(ln, x<-"ln")
match_mp <- match(mp, x<-"mp")
match_mn <- match(mn, x<-"mn")

compound_comp_exradii <- bind_rows(match_lp, match_ln, match_mp, match_mn)

all <- combine_tables(lp, ln, mp, mn)

filtered_exradii_v_rad <- filtering(all)

filtered_exradii_v_rad <- filter(filtered_exradii_v_rad, overlap > 0 & pathway_size >0)

# Creating rad v drug changes -------------------------------------------------------------------------
all_comp <- bind_rows(filtered_rad_v_sham, filtered_veh_v_rad, filtered_exradi_v_rad, filtered_exradii_v_rad)
all_comp <- all_comp %>% 
  ungroup() %>%
  select(pathway) %>%
  distinct(pathway) %>%
  left_join(filtered_rad_v_sham, by = "pathway") %>%
  left_join(filtered_veh_v_rad, by = "pathway") %>%
  left_join(filtered_exradi_v_rad, by = "pathway") %>%
  left_join(filtered_exradii_v_rad, by = "pathway") %>%
  select(pathway, overlap.x, pathway_size.x, p_value.x, overlap.y, pathway_size.y, p_value.y, overlap.x.x, pathway_size.x.x, p_value.x.x, overlap.y.y, pathway_size.y.y, p_value.y.y) %>%
  mutate(rad_v_sham = paste(overlap.x, pathway_size.x, " ", sep = '_'),
         veh_v_rad = paste(overlap.y, pathway_size.y, " ", sep = "_"), 
         exradi_v_rad = paste(overlap.x.x, pathway_size.x.x, " ", sep = "_"), 
         exradii_v_rad = paste(overlap.y.y, pathway_size.y.y, " ", sep = "_")) %>%
  select(pathway, rad_v_sham, p_value.x, veh_v_rad, p_value.y, exradi_v_rad, p_value.x.x, exradii_v_rad, p_value.y.y)

# Export tables ----------------------------------------------------------------------------------

write_tsv(match_ms2_pathway(compound_comp_rad_v_sham), paste0(project_dir, 
                                           "/data_tables/compound_comp_rad_v_sham.tsv"))

write_tsv(match_ms2_pathway(compound_comp_vehicle_v_rad), paste0(project_dir, 
                                       "/data_tables/compound_comp_vehicle_v_rad.tsv"))
write_tsv(match_ms2_pathway(compound_comp_exradi), paste0(project_dir, 
                                           "/data_tables/compound_comp_exradi.tsv"))
write_tsv(match_ms2_pathway(compound_comp_exradii), paste0(project_dir, 
                                           "/data_tables/compound_comp_exradii.tsv"))

write_tsv(all_comp, paste0(project_dir, "/data_tables/pathways_all.tsv"))
