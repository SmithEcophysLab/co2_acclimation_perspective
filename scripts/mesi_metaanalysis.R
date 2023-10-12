# MESI data trends for leaf/whole-plant responses to eCO2.
# Note: code assumes that the path to this script is the root working directory

## Libraries
library(tidyverse)
library(metafor)
library(MAd)

## Load MESI dataset (repository cloned from 
## https://github.com/MESI-organization/mesi-db)
df <- read.csv(
  #"../../mesi-db/data/mesi_main.csv"
  ) 

## See trait options
unique(df$response)
vars_of_interest <- c("anet", "jmax", "leaf_n_mass", 
                      "root_shoot_ratio", "total_biomass", "vcmax",
                      "soil_no3-n", "soil_nh4-n", "soil_in", "gs",
                      "root_n_uptake", "root_nh4_uptake", "root_no3_uptake",
                      "wue", "lai", "lai_max", "gpp", "asat")

## Assess treatment combinations. Searching for CO2 experiments only
unique(df$treatment)

## Filter MESI to include CO2-only experiments with traits of interest
df_filtered <- df %>%
  
  # Filter to include CO2-only experiments
  filter(treatment == "c") %>%
  
  # Remove experiments without defined medium
  #filter(!is.na(experiment_type)) %>%
  
  # Only include relevant measurements
  filter(response %in% vars_of_interest) %>%
  
  # Keep only essential columns
  dplyr::select(id, exp, response, treatment, sampling_year, 
                c_c, c_t, x_c, x_t, sd_c, sd_t, rep_c, rep_t)

## Calculate log response ratio and its variance, remove any lines with NA
log_resp <- df_filtered %>%
  mutate(across(c_c:rep_t, as.numeric)) %>%
  metafor::escalc( 
    measure = "ROM", 
    m1i = x_t, sd1i = sd_t, n1i = rep_t, 
    m2i = x_c, sd2i = sd_c, n2i = rep_c, 
    data = ., 
    append = TRUE, 
    var.names = c("logr", "logr_var")) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_t)) %>%
  filter(!is.na(logr_var) & !is.na(logr))

## Aggregate dependent effect sizes by experiment and response variable for
## meta-analysis
log_resp_agg <- log_resp %>%
  dplyr::select(-id) %>%
  mutate(id = paste(exp, response, sep = "_XXX_")) %>%
  agg(id = id, es = logr, var = logr_var, 
      cor = 1.0, method = "BHHR", data = .) %>%
  mutate(id = str_split(id, "_XXX_"),
         exp = purrr::map_chr(id, 1),
         response = purrr::map_chr(id, 2)) %>%
  dplyr::select(exp, response, logr = es, logr_var = var) %>%
  left_join(log_resp %>% group_by(exp, response, treatment) %>%
              summarize(rep_c = sum(rep_c), rep_t = sum(rep_t))) %>%
  mutate(logr_se = sqrt(logr_var) / sqrt(rep_c)) %>%
  mutate(response = ifelse(response %in% c("soil_no3-n", "soil_nh4-n",
                                           "soil_in"),
                           "available_n", response),
         response = ifelse(response %in% c("root_n_uptake", "root_nh4_uptake", 
                                           "root_no3_uptake"),
                           "n_uptake", response),
         response = ifelse(response %in% c("lai", "lai_max"), "lai", response)) %>%
  filter(logr > -2)

## Create new vars of interest vector with new "available_n", "n_uptake", and
## "lai" cols
vars_of_interest2 <- c("anet", "jmax", "leaf_n_mass", "gs",
                       "root_shoot_ratio", "total_biomass", "vcmax",
                       "available_n", "n_uptake", "lai", "gpp", "asat")

## Pull analyse_meta.R script from Beni's repository to replicate approach used
## in Nwg meta-analysis

#source("../../lt_cn_review/R/analyse_meta.R")
out  <- purrr::map(as.list(vars_of_interest2), 
                   ~analyse_meta(log_resp_agg %>% 
                                   rename(var = response),
                                 nam_target = .))
names(out) <- vars_of_interest2
df_box <- purrr::map_dfr(out, "df_box") %>%
  mutate(var = factor(var, 
                      levels = c("available_n", "n_uptake", "root_shoot_ratio",
                                 "total_biomass", "lai", "gpp",
                                 "wue", "gs", "leaf_n_mass", "jmax", "vcmax", "asat", 
                                 "anet")))

## Subset only needed vars, factor to desirable plot order that 
## corresponds with leaf-whole plant-ecosystem levels
log_resp_agg_plot_subset <- log_resp_agg %>%
  filter(response %in% c("anet", "jmax", "leaf_n_mass", 
         "root_shoot_ratio", "total_biomass", "vcmax", "gs",
         "available_n", "lai", "gpp", "asat", "n_uptake"))

log_resp_agg_plot_subset$response <- factor(log_resp_agg_plot_subset$response,
                                levels = c("available_n", "n_uptake", 
                                           "root_shoot_ratio",
                                           "total_biomass", "lai", "gpp",
                                           "wue", "gs", "leaf_n_mass", "jmax", 
                                           "vcmax", "asat", 
                                           "anet"))

# Make plot
#png("../figs/co2_perspective_general_trends.png", height = 4.5, width = 8, 
#    units = "in", res = 600)                               
ggplot(data = log_resp_agg_plot_subset, aes(x = response, y = logr)) +
  geom_jitter(aes(size = logr_se), alpha = 0.4, width = 0.20, shape = 21) +
  geom_errorbar(data = df_box, 
                aes(x = var, y = middle, ymin = ymin, ymax = ymax), 
                width = 0.5, size = 1) +
  geom_point(data = df_box, aes(x = var, y = middle), 
             size = 2, shape = 21, fill = "black") +
  # geom_crossbar(data = df_box, 
  #               aes(x = var, y = middle, ymin = ymin, ymax = ymax),
  #               fill = "turquoise", color = "black", alpha = 0.6, width = 0.5) +
  scale_x_discrete(labels = c("Inorganic N availability", 
                              "Root N uptake",
                              "Root:shoot", 
                              "Total biomass", 
                              "LAI", 
                              "GPP",
                              expression("g"["s"]),
                              expression("N"["mass"]),
                              expression("J"["max"]), 
                              expression("V"["cmax"]),
                              expression("A"["sat"]),
                              expression("A"["net"]))) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "red") +
  labs(x = "Response variable", y = expression("Response to elevated CO"["2"]),
       size = "Standard error of\nlog response ratio") + 
  scale_y_continuous(limits = c(-0.8, 1.2), 
                     breaks = c(-0.8, -0.4, 0.0, 0.4, 0.8, 1.2)) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major.y = element_blank())
#dev.off()







