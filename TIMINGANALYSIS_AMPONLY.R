# Timing of copy-number gains relative to clock-like point mutations (Fig. 4c).
# Sourced by 11_fig4.R: uses DATA_DIR and the clock-like rate (mutation.lm.all) from that session.
# Molecular time (years) = early SNV rate per Mb x 3,140 Mb (genome size) / clock-like mutations per year.

setwd(file.path(DATA_DIR, '50_additional/timing_analysis/timing_output'))

library(dplyr)
library(data.table)
library(purrr)
library(stringr)

# regression coefficient for age at diagnosis (assumed to exist in the session)
diag_coef <- as.numeric(mutation.lm.all$coefficients[["Diagnosis_age"]])

# function for one patient
run_one_patient <- function(patient, root = getwd()) {
  files <- list.files(path = root,
                      pattern = paste0("^", patient, ".*earlysnv$"),
                      full.names = TRUE) 
  if (patient == 'C10') {files <- files[!grepl('META', files)]}
  if (length(files) == 0) return(NULL)
  
  # read files -> one table
  df <- map(files, function(x) {
    read.table(x,
               col.names = c("broad_amp", "dist", "early_SNV", "early_SNV_range",
                             "early_SNVrate/Mb", "early_SNVrate_range", "early_SNV_var"),
               check.names = TRUE,  # "early_SNVrate/Mb" -> "early_SNVrate.Mb"
               stringsAsFactors = FALSE) %>%
      mutate(sample  = sub("_clean.*", "", basename(x)),
             patient = patient)
  }) %>% rbindlist(fill = TRUE) %>% as_tibble()
  
  # note: with check.names = TRUE the column is available as early_SNVrate.Mb
  
  # annotation
  df_annot <- df %>%
    group_by(patient) %>%
    mutate(total_samples = n_distinct(sample)) %>%
    group_by(patient, broad_amp) %>%
    mutate(
      n_samples   = n_distinct(sample),
      samples_str = paste(sort(unique(sample)), collapse = ";"),
      mean_rate   = mean(early_SNVrate.Mb, na.rm = TRUE),
      sd_rate     = sd(early_SNVrate.Mb,   na.rm = TRUE),
      molecular_time_mean = mean_rate * 3140 / diag_coef,
      molecular_time_sd   = sd_rate   * 3140 / diag_coef,
      event       = if_else(n_samples == total_samples, "MRCA", "branch")
    ) %>%
    ungroup() %>%
    # important: keep one row per patient and amplification (not shared across patients)
    distinct(patient, broad_amp, .keep_all = TRUE) %>%
    # drop raw columns that are not needed
    select(-sample, -early_SNV, -early_SNV_range, -early_SNVrate.Mb,
           -early_SNV_var, -early_SNVrate_range, -dist, -mean_rate, -sd_rate)
  
  df_annot
}

# run for several patients and combine
patients <- c("C04","C01","C09","C03","C08","C10")

# run + combine
all_annot <- map(patients, run_one_patient) %>%
  bind_rows() %>%
  arrange(patient, broad_amp) %>% 
  filter(n_samples > 1)

# check the result
all_annot

### PRUNING BASED ON THE COPY NUMBER PHYLOGENY 
all_annot <- all_annot |> 
  filter(!(broad_amp == '17pAMP' & patient == 'C01'),
         !(broad_amp %in% c('21pAMP', '9pAMP') & patient == 'C04'),
         !(broad_amp != '14qAMP' & patient == 'C08'))

# helper that builds a per-patient summary
summarize_patient <- function(df, patient_id, event_name, exclude_amps = character()) {
  sub <- df %>% filter(patient == patient_id, !broad_amp %in% exclude_amps)
  if (nrow(sub) == 0) return(NULL)
  
  # merge samples_str (unique, joined with ';')
  all_samples <- sort(unique(unlist(strsplit(paste(sub$samples_str, collapse = ";"), ";"))))
  
  tibble(
    broad_amp           = "WGD",
    patient             = patient_id,
    total_samples       = dplyr::first(sub$total_samples),
    n_samples           = length(all_samples),
    samples_str         = paste(all_samples, collapse = ";"),
    molecular_time_mean = mean(sub$molecular_time_mean, na.rm = TRUE),
    molecular_time_sd   = sd  (sub$molecular_time_mean, na.rm = TRUE),  # sd of the means
    event               = event_name
  )
}

# summary rows
sum_c01 <- summarize_patient(all_annot,event_name = 'branch', "C01")                 # all included
sum_c09 <- summarize_patient(all_annot,event_name = 'MRCA', "C09", "21WCA")        # 21WCA excluded
sum_c10 <- summarize_patient(all_annot,event_name = 'MRCA', "C10")                 # all included

# remaining patients (C01, C09, C10 removed); the 21WCA exclusion for C09 applies only to the summary above
rest <- all_annot %>%
  filter(!patient %in% c("C01", "C09", "C10"))

# combine
final_annot <- bind_rows(sum_c01, sum_c09, sum_c10, rest) %>%
  arrange(patient, broad_amp)

final_annot
