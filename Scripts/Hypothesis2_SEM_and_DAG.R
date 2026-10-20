library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(brms)
library(cmdstanr)
library(posterior)
library(bayestestR)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


exp_data<-read.csv("metacom_exp_data_final.csv")

exp_data$metacommunity <- factor(exp_data$metacommunity)
exp_data$wattage <- factor(exp_data$wattage, levels=c("0","100","200","300"))
exp_data$dispersal <- factor(exp_data$dispersal, levels=c("none","low","high"))
exp_data$tank <- factor(exp_data$tank)
exp_data$sample_date <- as.Date(exp_data$sample_date, format = "%m/%d/%y") #sometimes this is necessary
exp_data$sample_date <- as.character(exp_data$sample_date)


# subset and reshape data
df_sem_subset0 <- exp_data %>%
  dplyr::select(
    ln_GPP_umol_L_h,CTPI_scaled,slope_size_spectra,
    ln_fluoro_ug_C_L, inv_kTc,
    totalN,
    dispersal, metacommunity, tank, ln_zoop_ug_L, week_number
  ) %>%
  tidyr::drop_na()

df_sem_subset <- df_sem_subset0 %>%
  mutate(
    inv_kTc             = inv_kTc,
    size_spec        = slope_size_spectra,
    thermaltraits        = CTPI_scaled, # = Community Thermal Phenotype Index
    ln_fluoro_ug_C_L  = ln_fluoro_ug_C_L,
    ln_totalN          = scale(log(totalN)),
    ln_GPP_umol_L_h   = as.numeric(ln_GPP_umol_L_h),
    dispersal           = factor(dispersal),
    metacommunity       = factor(metacommunity),
    tank                = factor(tank),
    zoop_ug_C_L       = scale(ln_zoop_ug_L),
    week_number         = factor(week_number, levels = c("4","8","12"), ordered = FALSE)
  )

# Each submodel now includes week interactions on all fixed effects
# This allows every path slope to vary by week

bf_size_wk <- bf(
  size_spec ~ (dispersal + inv_kTc + zoop_ug_C_L) * week_number +
    (1 | metacommunity),
  family = gaussian()
)

bf_comp_wk <- bf(
  thermaltraits ~ (dispersal + inv_kTc + zoop_ug_C_L) * week_number +
    (1 | metacommunity),
  family = gaussian()
)

bf_biomass_wk <- bf(
  ln_fluoro_ug_C_L ~ (dispersal + inv_kTc + ln_totalN + zoop_ug_C_L) * week_number +
    (1 | metacommunity)
)

bf_zoop_wk <- bf(
  zoop_ug_C_L ~ (inv_kTc + dispersal) * week_number +
    (1 | metacommunity),
  family = gaussian()
)

bf_gpp_wk <- bf(
  ln_GPP_umol_L_h ~ (size_spec + thermaltraits + ln_fluoro_ug_C_L +
                       inv_kTc + ln_totalN) * week_number +
    (1 | metacommunity),
  family = gaussian()
)

# set priors (mostly uninformative except temperature and mass terms from metabolic theory)
priors_all <- c(
  # Intercepts
  set_prior("normal(0,1)", class = "Intercept", resp = "sizespec"),   
  set_prior("normal(0,1)", class = "Intercept", resp = "thermaltraits"),
  set_prior("normal(0,1)", class = "Intercept", resp = "lnfluorougCL"),
  set_prior("normal(0,1)", class = "Intercept", resp = "lnGPPumolLh"),
  set_prior("normal(0,1)", class = "Intercept", resp = "zoopugCL"),
  
  # Phytoplankton size spectra
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc",       resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow",  resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh", resp = "sizespec"),
  set_prior("normal(0.5,0.5)", class = "b", coef = "zoop_ug_C_L", resp = "sizespec"),
  
  # Phytoplankton thermal trait composition (Community Thermal Phenotype Index)
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc",       resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow",  resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh", resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "zoop_ug_C_L",   resp = "thermaltraits"),
  
  # Phyto biomass (fluoroprobe-measured grams of Carbon per Litre)
  set_prior("normal(-0.5,0.5)", class = "b", coef = "inv_kTc",       resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)",    class = "b", coef = "dispersallow",  resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)",    class = "b", coef = "dispersalhigh", resp = "lnfluorougCL"),
  set_prior("normal(1,0.5)",    class = "b", coef = "ln_totalN",     resp = "lnfluorougCL"),
  set_prior("normal(-0.5,0.5)", class = "b", coef = "zoop_ug_C_L",   resp = "lnfluorougCL"),
  
  # Zoop biomass
  set_prior("normal(0.5,0.5)", class = "b", coef = "inv_kTc",       resp = "zoopugCL"),
  set_prior("normal(0,0.5)",   class = "b", coef = "dispersallow",  resp = "zoopugCL"),
  set_prior("normal(0,0.5)",   class = "b", coef = "dispersalhigh", resp = "zoopugCL"),
  
  # GPP 
  set_prior("normal(0.65,0.28)", class = "b", coef = "inv_kTc",          resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)",     class = "b", coef = "ln_totalN",        resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)",     class = "b", coef = "size_spec",        resp = "lnGPPumolLh"), 
  set_prior("normal(0,0.5)",     class = "b", coef = "thermaltraits",        resp = "lnGPPumolLh"), 
  set_prior("normal(0.75,0.5)",  class = "b", coef = "ln_fluoro_ug_C_L", resp = "lnGPPumolLh"),
  
  # Residual SDs 
  set_prior("exponential(1)", class = "sigma", resp = "sizespec"),
  set_prior("exponential(1)", class = "sigma", resp = "thermaltraits"),
  set_prior("exponential(1)", class = "sigma", resp = "lnfluorougCL"),
  set_prior("exponential(1)", class = "sigma", resp = "lnGPPumolLh"),
  set_prior("exponential(1)", class = "sigma", resp = "zoopugCL")
)


priors_wk <- c(
  priors_all,
  set_prior("normal(0,0.5)", class = "b", coef = "week_number8",  resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "week_number12", resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "week_number8",  resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "week_number12", resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "week_number8",  resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "week_number12", resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "week_number8",  resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "week_number12", resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "week_number8",  resp = "zoopugCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "week_number12", resp = "zoopugCL"),
  
  # Interaction terms: how each path slope changes from week 4 to 8 or 12
  # GPP submodel interactions
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number8",          resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number12",         resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "size_spec:week_number8",        resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "size_spec:week_number12",       resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "thermaltraits:week_number8",    resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "thermaltraits:week_number12",   resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "ln_fluoro_ug_C_L:week_number8", resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "ln_fluoro_ug_C_L:week_number12",resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "ln_totalN:week_number8",        resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)", class = "b", coef = "ln_totalN:week_number12",       resp = "lnGPPumolLh"),
  
  # Biomass submodel interactions
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number8",          resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number12",         resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow:week_number8",     resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow:week_number12",    resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh:week_number8",    resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh:week_number12",   resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "ln_totalN:week_number8",        resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "ln_totalN:week_number12",       resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "zoop_ug_C_L:week_number8",      resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "zoop_ug_C_L:week_number12",     resp = "lnfluorougCL"),
  
  # Thermal traits submodel interactions
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number8",          resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number12",         resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow:week_number8",     resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow:week_number12",    resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh:week_number8",    resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh:week_number12",   resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "zoop_ug_C_L:week_number8",      resp = "thermaltraits"),
  set_prior("normal(0,0.5)", class = "b", coef = "zoop_ug_C_L:week_number12",     resp = "thermaltraits"),
  
  # Size spectra submodel interactions
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number8",          resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number12",         resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow:week_number8",     resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow:week_number12",    resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh:week_number8",    resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh:week_number12",   resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "zoop_ug_C_L:week_number8",      resp = "sizespec"),
  set_prior("normal(0,0.5)", class = "b", coef = "zoop_ug_C_L:week_number12",     resp = "sizespec"),
  
  # Zoop submodel interactions
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number8",          resp = "zoopugCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc:week_number12",         resp = "zoopugCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow:week_number8",     resp = "zoopugCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow:week_number12",    resp = "zoopugCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh:week_number8",    resp = "zoopugCL"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh:week_number12",   resp = "zoopugCL")
)

# Fit the model
fit_wk <- brm(
  bf_size_wk + bf_comp_wk + bf_biomass_wk + bf_gpp_wk + bf_zoop_wk +
    set_rescor(FALSE),
  data      = df_sem_subset,
  prior     = priors_wk,
  chains    = 4, cores = 4, iter = 5000,
  backend   = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  control   = list(adapt_delta = 0.99, max_treedepth = 15),
  seed      = 2025
)

# Exctract week-specific path estimates
draws_wk <- posterior::as_draws_df(fit_wk)
get_wk <- function(var) as.numeric(draws_wk[[var]])

# Helper: week-specific slope = main effect + week interaction
path_at_week <- function(main_coef, int_coef_8, int_coef_12) {
  main  <- get_wk(main_coef)
  d8    <- if (!is.null(int_coef_8)  && int_coef_8  %in% names(draws_wk))
    get_wk(int_coef_8)  else rep(0, nrow(draws_wk))
  d12   <- if (!is.null(int_coef_12) && int_coef_12 %in% names(draws_wk))
    get_wk(int_coef_12) else rep(0, nrow(draws_wk))
  list(wk4 = main,
       wk8 = main + d8,
       wk12 = main + d12)
}

# Direct paths to GPP 
temp_gpp   <- path_at_week("b_lnGPPumolLh_inv_kTc",
                           "b_lnGPPumolLh_inv_kTc:week_number8",
                           "b_lnGPPumolLh_inv_kTc:week_number12")
size_gpp   <- path_at_week("b_lnGPPumolLh_size_spec",
                           "b_lnGPPumolLh_size_spec:week_number8",
                           "b_lnGPPumolLh_size_spec:week_number12")
comp_gpp   <- path_at_week("b_lnGPPumolLh_thermaltraits",
                           "b_lnGPPumolLh_thermaltraits:week_number8",
                           "b_lnGPPumolLh_thermaltraits:week_number12")
phyto_gpp  <- path_at_week("b_lnGPPumolLh_ln_fluoro_ug_C_L",
                           "b_lnGPPumolLh_ln_fluoro_ug_C_L:week_number8",
                           "b_lnGPPumolLh_ln_fluoro_ug_C_L:week_number12")
N_gpp      <- path_at_week("b_lnGPPumolLh_ln_totalN",
                           "b_lnGPPumolLh_ln_totalN:week_number8",
                           "b_lnGPPumolLh_ln_totalN:week_number12")

# Paths into mediators of GPP
temp_phyto <- path_at_week("b_lnfluorougCL_inv_kTc",
                           "b_lnfluorougCL_inv_kTc:week_number8",
                           "b_lnfluorougCL_inv_kTc:week_number12")
temp_comp  <- path_at_week("b_thermaltraits_inv_kTc",
                           "b_thermaltraits_inv_kTc:week_number8",
                           "b_thermaltraits_inv_kTc:week_number12")
temp_size  <- path_at_week("b_sizespec_inv_kTc",
                           "b_sizespec_inv_kTc:week_number8",
                           "b_sizespec_inv_kTc:week_number12")
temp_zoop  <- path_at_week("b_zoopugCL_inv_kTc",
                           "b_zoopugCL_inv_kTc:week_number8",
                           "b_zoopugCL_inv_kTc:week_number12")

N_phyto    <- path_at_week("b_lnfluorougCL_ln_totalN",
                           "b_lnfluorougCL_ln_totalN:week_number8",
                           "b_lnfluorougCL_ln_totalN:week_number12")

high_phyto <- path_at_week("b_lnfluorougCL_dispersalhigh",
                           "b_lnfluorougCL_dispersalhigh:week_number8",
                           "b_lnfluorougCL_dispersalhigh:week_number12")
low_phyto  <- path_at_week("b_lnfluorougCL_dispersallow",
                           "b_lnfluorougCL_dispersallow:week_number8",
                           "b_lnfluorougCL_dispersallow:week_number12")
high_comp  <- path_at_week("b_thermaltraits_dispersalhigh",
                           "b_thermaltraits_dispersalhigh:week_number8",
                           "b_thermaltraits_dispersalhigh:week_number12")
low_comp   <- path_at_week("b_thermaltraits_dispersallow",
                           "b_thermaltraits_dispersallow:week_number8",
                           "b_thermaltraits_dispersallow:week_number12")
high_size  <- path_at_week("b_sizespec_dispersalhigh",
                           "b_sizespec_dispersalhigh:week_number8",
                           "b_sizespec_dispersalhigh:week_number12")
low_size   <- path_at_week("b_sizespec_dispersallow",
                           "b_sizespec_dispersallow:week_number8",
                           "b_sizespec_dispersallow:week_number12")
high_zoop  <- path_at_week("b_zoopugCL_dispersalhigh",
                           "b_zoopugCL_dispersalhigh:week_number8",
                           "b_zoopugCL_dispersalhigh:week_number12")
low_zoop   <- path_at_week("b_zoopugCL_dispersallow",
                           "b_zoopugCL_dispersallow:week_number8",
                           "b_zoopugCL_dispersallow:week_number12")

zoop_size  <- path_at_week("b_sizespec_zoop_ug_C_L",
                           "b_sizespec_zoop_ug_C_L:week_number8",
                           "b_sizespec_zoop_ug_C_L:week_number12")
zoop_comp  <- path_at_week("b_thermaltraits_zoop_ug_C_L",
                           "b_thermaltraits_zoop_ug_C_L:week_number8",
                           "b_thermaltraits_zoop_ug_C_L:week_number12")
zoop_phyto <- path_at_week("b_lnfluorougCL_zoop_ug_C_L",
                           "b_lnfluorougCL_zoop_ug_C_L:week_number8",
                           "b_lnfluorougCL_zoop_ug_C_L:week_number12")

# Indirect effects (mediators by week)
for (wk in c("wk4", "wk8", "wk12")) {
  assign(paste0("ind_high_", wk),
         high_size[[wk]] * size_gpp[[wk]] +
           high_comp[[wk]] * comp_gpp[[wk]] +
           high_phyto[[wk]] * phyto_gpp[[wk]])
  
  assign(paste0("ind_low_", wk),
         low_size[[wk]] * size_gpp[[wk]] +
           low_comp[[wk]] * comp_gpp[[wk]] +
           low_phyto[[wk]] * phyto_gpp[[wk]])
  
  assign(paste0("ind_temp_", wk),
         temp_size[[wk]] * size_gpp[[wk]] +
           temp_comp[[wk]] * comp_gpp[[wk]] +
           temp_phyto[[wk]] * phyto_gpp[[wk]])
}

# Set ROPE thresholds per week, 20% SD(GPP)

thr_wk4  <- 0.20 * sd(df_sem_subset %>% filter(week_number == "4")  %>% pull(ln_GPP_umol_L_h), na.rm = TRUE)
thr_wk8  <- 0.20 * sd(df_sem_subset %>% filter(week_number == "8")  %>% pull(ln_GPP_umol_L_h), na.rm = TRUE)
thr_wk12 <- 0.20 * sd(df_sem_subset %>% filter(week_number == "12") %>% pull(ln_GPP_umol_L_h), na.rm = TRUE)

# set pooled ROPE threshold: 20% SD(GPP) across all weeks. used for direct/indirect effect table to ensure fair comparison across weeks (because SD(GPP) differs across weeks)

thr_pooled <- 0.20 * sd(df_sem_subset %>% pull(ln_GPP_umol_L_h), na.rm = TRUE)

# SUMMARY TABLE: direct and indirect effects per week and "significance" based on % of 89% HDI within ROPE

make_row <- function(draws, thr, path, week) {
  tibble(
    Path       = path,
    Week       = week,
    Coefficient = round(mean(draws), 3),
    L89        = round(quantile(draws, 0.055), 3),
    U89        = round(quantile(draws, 0.945), 3),
    HDIinROPE  = round(hdi_in_rope_pct(draws, thr), 1),
    Decision   = case_when(
      hdi_in_rope_pct(draws, thr) <= 5  ~ "nonzero",
      hdi_in_rope_pct(draws, thr) >= 95 ~ "null",
      TRUE ~ "ambiguous"
    )
  )
}

direct_indirect_wk <- bind_rows(
  make_row(ind_high_wk4,  thr_pooled, "High dispersal → GPP (indirect)", "4"),
  make_row(ind_high_wk8,  thr_pooled, "High dispersal → GPP (indirect)", "8"),
  make_row(ind_high_wk12, thr_pooled, "High dispersal → GPP (indirect)", "12"),
  make_row(ind_low_wk4,   thr_pooled, "Low dispersal → GPP (indirect)",  "4"),
  make_row(ind_low_wk8,   thr_pooled, "Low dispersal → GPP (indirect)",  "8"),
  make_row(ind_low_wk12,  thr_pooled, "Low dispersal → GPP (indirect)",  "12"),
  make_row(temp_gpp$wk4,  thr_pooled, "Temp → GPP (direct)",             "4"),
  make_row(temp_gpp$wk8,  thr_pooled, "Temp → GPP (direct)",             "8"),
  make_row(temp_gpp$wk12, thr_pooled, "Temp → GPP (direct)",             "12"),
  make_row(ind_temp_wk4,  thr_pooled, "Temp → GPP (indirect)",           "4"),
  make_row(ind_temp_wk8,  thr_pooled, "Temp → GPP (indirect)",           "8"),
  make_row(ind_temp_wk12, thr_pooled, "Temp → GPP (indirect)",           "12")
)
### RESULTS IN TABLE 1:
print(direct_indirect_wk, n = Inf)

##### ROPE sensitivity version:

sd_pooled <- sd(df_sem_subset %>% pull(ln_GPP_umol_L_h), na.rm = TRUE)

make_row_sens <- function(draws, path, week) {
  widths <- c(0.05, 0.10, 0.20)
  pcts   <- sapply(widths, function(w) {
    round(hdi_in_rope_pct(draws, w * sd_pooled, ci = 0.89), 1)
  })
  tibble(
    Path        = path,
    Week        = week,
    Coefficient = round(mean(draws), 3),
    L89         = round(quantile(draws, 0.055), 3),
    U89         = round(quantile(draws, 0.945), 3),
    ROPE_5_pct  = pcts[1],
    ROPE_10_pct = pcts[2],
    ROPE_20_pct = pcts[3],
    sensitive   = (pcts[1] <= 5 & (pcts[2] > 5 | pcts[3] > 5)) |
      (pcts[2] <= 5 &  pcts[3] > 5)
  )
}

direct_indirect_sens <- bind_rows(
  make_row_sens(ind_high_wk4,  "High dispersal → GPP (indirect)", "4"),
  make_row_sens(ind_high_wk8,  "High dispersal → GPP (indirect)", "8"),
  make_row_sens(ind_high_wk12, "High dispersal → GPP (indirect)", "12"),
  make_row_sens(ind_low_wk4,   "Low dispersal → GPP (indirect)",  "4"),
  make_row_sens(ind_low_wk8,   "Low dispersal → GPP (indirect)",  "8"),
  make_row_sens(ind_low_wk12,  "Low dispersal → GPP (indirect)",  "12"),
  make_row_sens(temp_gpp$wk4,  "Temp → GPP (direct)",             "4"),
  make_row_sens(temp_gpp$wk8,  "Temp → GPP (direct)",             "8"),
  make_row_sens(temp_gpp$wk12, "Temp → GPP (direct)",             "12"),
  make_row_sens(ind_temp_wk4,  "Temp → GPP (indirect)",           "4"),
  make_row_sens(ind_temp_wk8,  "Temp → GPP (indirect)",           "8"),
  make_row_sens(ind_temp_wk12, "Temp → GPP (indirect)",           "12")
)

# Table of direct and indirect sensitivity to ROPE width for overall effects on GPP only. 
# Shows sensitivity to ROPE width (5/10/20% SD(GPP)) of paths presented in main text Table 1

print(direct_indirect_sens, n = Inf)

# TABLE OF ALL EDGE POSTERIORS (Table S3 in Supporting Info)

# function building tribble by week
build_edge_summ <- function(wk_label, thr_gpp, thr_size, thr_comp,
                            thr_phyto, thr_zoop) {
  wk <- paste0("wk", wk_label)
  
  edges <- tibble::tribble(
    ~from,               ~to,              ~draws,
    "Temp",              "GPP",            temp_gpp[[wk]],
    "Total N",           "GPP",            N_gpp[[wk]],
    "Size spectra",      "GPP",            size_gpp[[wk]],
    "Thermal trait comp","GPP",            comp_gpp[[wk]],
    "Phyto biomass",     "GPP",            phyto_gpp[[wk]],
    "Temp",              "Zoop biomass",   temp_zoop[[wk]],
    "Temp",              "Size spectra",   temp_size[[wk]],
    "Temp",              "Thermal trait comp", temp_comp[[wk]],
    "Temp",              "Phyto biomass",  temp_phyto[[wk]],
    "Total N",           "Phyto biomass",  N_phyto[[wk]],
    "Dispersal: low",    "Zoop biomass",   low_zoop[[wk]],
    "Dispersal: high",   "Zoop biomass",   high_zoop[[wk]],
    "Dispersal: low",    "Size spectra",   low_size[[wk]],
    "Dispersal: high",   "Size spectra",   high_size[[wk]],
    "Dispersal: low",    "Thermal trait comp", low_comp[[wk]],
    "Dispersal: high",   "Thermal trait comp", high_comp[[wk]],
    "Dispersal: low",    "Phyto biomass",  low_phyto[[wk]],
    "Dispersal: high",   "Phyto biomass",  high_phyto[[wk]],
    "Zoop biomass",      "Size spectra",   zoop_size[[wk]],
    "Zoop biomass",      "Thermal trait comp", zoop_comp[[wk]],
    "Zoop biomass",      "Phyto biomass",  zoop_phyto[[wk]]
  ) %>%
    mutate(
      week      = wk_label,
      Mean      = map_dbl(draws, mean),
      L89       = map_dbl(draws, ~quantile(.x, 0.055)),
      U89       = map_dbl(draws, ~quantile(.x, 0.945)),
      thr       = case_when(
        to == "GPP"                ~ thr_gpp,
        to == "Size spectra"       ~ thr_size,
        to == "Thermal trait comp" ~ thr_comp,
        to == "Phyto biomass"      ~ thr_phyto,
        to == "Zoop biomass"       ~ thr_zoop,
        TRUE                       ~ NA_real_
      ),
      HDIinROPE = map2_dbl(draws, thr, ~hdi_in_rope_pct(.x, .y, ci = 0.89)),
      decision  = case_when(
        HDIinROPE >= 95 ~ "null",
        HDIinROPE <=  5 ~ "nonzero",
        TRUE            ~ "ambiguous"
      )
    )
  edges
}

# ROPE thresholds per mediator per week 
sd_size4  <- sd(df_sem_subset %>% filter(week_number=="4")  %>% pull(size_spec),         na.rm=TRUE)
sd_comp4  <- sd(df_sem_subset %>% filter(week_number=="4")  %>% pull(thermaltraits),     na.rm=TRUE)
sd_phyto4 <- sd(df_sem_subset %>% filter(week_number=="4")  %>% pull(ln_fluoro_ug_C_L), na.rm=TRUE)
sd_zoop4  <- sd(df_sem_subset %>% filter(week_number=="4")  %>% pull(zoop_ug_C_L),      na.rm=TRUE)
sd_size8  <- sd(df_sem_subset %>% filter(week_number=="8")  %>% pull(size_spec),         na.rm=TRUE)
sd_comp8  <- sd(df_sem_subset %>% filter(week_number=="8")  %>% pull(thermaltraits),     na.rm=TRUE)
sd_phyto8 <- sd(df_sem_subset %>% filter(week_number=="8")  %>% pull(ln_fluoro_ug_C_L), na.rm=TRUE)
sd_zoop8  <- sd(df_sem_subset %>% filter(week_number=="8")  %>% pull(zoop_ug_C_L),      na.rm=TRUE)
sd_size12 <- sd(df_sem_subset %>% filter(week_number=="12") %>% pull(size_spec),         na.rm=TRUE)
sd_comp12 <- sd(df_sem_subset %>% filter(week_number=="12") %>% pull(thermaltraits),     na.rm=TRUE)
sd_phyto12<- sd(df_sem_subset %>% filter(week_number=="12") %>% pull(ln_fluoro_ug_C_L), na.rm=TRUE)
sd_zoop12 <- sd(df_sem_subset %>% filter(week_number=="12") %>% pull(zoop_ug_C_L),      na.rm=TRUE)

thr_pooled <- 0.20 * sd(df_sem_subset %>% pull(ln_GPP_umol_L_h), na.rm = TRUE)

# Then rebuild edge summaries using thr_pooled for GPP paths in all weeks
edge_summ_wk4  <- build_edge_summ("4",  thr_pooled, 0.20*sd_size4,
                                  0.20*sd_comp4, 0.20*sd_phyto4, 0.20*sd_zoop4)
edge_summ_wk8  <- build_edge_summ("8",  thr_pooled, 0.20*sd_size8,
                                  0.20*sd_comp8, 0.20*sd_phyto8, 0.20*sd_zoop8)
edge_summ_wk12 <- build_edge_summ("12", thr_pooled, 0.20*sd_size12,
                                  0.20*sd_comp12, 0.20*sd_phyto12, 0.20*sd_zoop12)

edge_summ_all_wk <- bind_rows(edge_summ_wk4, edge_summ_wk8, edge_summ_wk12)

# Table S3, supporting info
edge_summ_all_wk %>%
  select(from, to, week, Mean, L89, U89, HDIinROPE, decision) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  arrange(from, to, week) %>%
  print(n = Inf)

# Week interaction terms:
# Do path slopes significantly change from week 4 to 8 or 12?

interaction_terms <- tibble::tribble(
  ~Path,                              ~Contrast,        ~draws,
  "Temp → GPP",                       "Week 8 - Week 4",  get_wk("b_lnGPPumolLh_inv_kTc:week_number8"),
  "Temp → GPP",                       "Week 12 - Week 4", get_wk("b_lnGPPumolLh_inv_kTc:week_number12"),
  "Phyto biomass → GPP",              "Week 8 - Week 4",  get_wk("b_lnGPPumolLh_ln_fluoro_ug_C_L:week_number8"),
  "Phyto biomass → GPP",              "Week 12 - Week 4", get_wk("b_lnGPPumolLh_ln_fluoro_ug_C_L:week_number12"),
  "Thermal traits → GPP",             "Week 8 - Week 4",  get_wk("b_lnGPPumolLh_thermaltraits:week_number8"),
  "Thermal traits → GPP",             "Week 12 - Week 4", get_wk("b_lnGPPumolLh_thermaltraits:week_number12"),
  "Disp: high → Phyto biomass",       "Week 8 - Week 4",  get_wk("b_lnfluorougCL_dispersalhigh:week_number8"),
  "Disp: high → Phyto biomass",       "Week 12 - Week 4", get_wk("b_lnfluorougCL_dispersalhigh:week_number12"),
  "Disp: high → Thermal trait comp",  "Week 8 - Week 4",  get_wk("b_thermaltraits_dispersalhigh:week_number8"),
  "Disp: high → Thermal trait comp",  "Week 12 - Week 4", get_wk("b_thermaltraits_dispersalhigh:week_number12"),
  "Temp → Phyto biomass",             "Week 8 - Week 4",  get_wk("b_lnfluorougCL_inv_kTc:week_number8"),
  "Temp → Phyto biomass",             "Week 12 - Week 4", get_wk("b_lnfluorougCL_inv_kTc:week_number12"),
  "Temp → Thermal trait comp",        "Week 8 - Week 4",  get_wk("b_thermaltraits_inv_kTc:week_number8"),
  "Temp → Thermal trait comp",        "Week 12 - Week 4", get_wk("b_thermaltraits_inv_kTc:week_number12")
) %>%
  mutate(
    Mean      = map_dbl(draws, mean),
    L89       = map_dbl(draws, ~quantile(.x, 0.055)),
    U89       = map_dbl(draws, ~quantile(.x, 0.945)),
    HDIinROPE = map_dbl(draws, ~hdi_in_rope_pct(.x, thr_wk4, ci = 0.89)),
    Decision  = case_when(
      HDIinROPE <= 5  ~ "nonzero",
      HDIinROPE >= 95 ~ "null",
      TRUE            ~ "ambiguous"
    )
  ) %>%
  select(Path, Contrast, Mean, L89, U89, HDIinROPE, Decision) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

# this table explicitly tests temporal changes in SEM path coefficients: 
# Nonzero = week 8/12 path slope differed significantly from week 4
print(interaction_terms, n = Inf)


# ============================================================
# ROPE SENSITIVITY TABLE from week-interaction model
# Uses the same rope_sensitivity_df function as original workflow
# Applied separately to each week's edge summary
# ============================================================
sd_pooled <- sd(df_sem_subset %>% pull(ln_GPP_umol_L_h), na.rm = TRUE)

rope_sensitivity_df_draws_pooled <- function(edge_df, sdy) {
  widths <- c(0.05, 0.10, 0.20)
  purrr::map_dfr(
    seq_len(nrow(edge_df)), \(i) {
      draws <- edge_df$draws[[i]]
      purrr::map_dfc(widths, \(w) {
        thr <- w * sdy          # sdy is now sd_pooled, same for all weeks
        pct <- hdi_in_rope_pct(draws, thr, ci = 0.89)
        tibble(!!paste0("ROPE_", w*100, "_pct") := pct)
      }) |>
        mutate(Effect = paste(edge_df$from[i], "\u2192", edge_df$to[i]))
    }
  )
}

rope_sens_wk4_pooled  <- rope_sensitivity_df_draws_pooled(edge_summ_wk4,  sd_pooled)
rope_sens_wk8_pooled  <- rope_sensitivity_df_draws_pooled(edge_summ_wk8,  sd_pooled)
rope_sens_wk12_pooled <- rope_sensitivity_df_draws_pooled(edge_summ_wk12, sd_pooled)

rope_sensitivity_wk_pooled <- bind_rows(
  rope_sens_wk4_pooled  %>% mutate(Week = "4"),
  rope_sens_wk8_pooled  %>% mutate(Week = "8"),
  rope_sens_wk12_pooled %>% mutate(Week = "12")
) %>%
  mutate(
    sensitive_to_rope = (ROPE_5_pct < 5 & (ROPE_10_pct >= 5 | ROPE_20_pct >= 5)) |
      (ROPE_10_pct < 5 & ROPE_20_pct >= 5)
  ) %>%
  select(Effect, Week, ROPE_5_pct, ROPE_10_pct, ROPE_20_pct, sensitive_to_rope) %>%
  arrange(Effect, Week)

# Full ROPE sensitivity table for all paths in SEMs by week, corresponds with Table S4 in Supporting Info
print(rope_sensitivity_wk_pooled, n = Inf)


# Plot Directed Acyclic Graphs (DAG) of SEM results as in Fig S9, Supporting Info

PW_CONSISTENT <- 8.0
PW_SOMETIME   <- 2.5
PW_NULL       <- 3.0
COL_SIG       <- "#222222"
COL_NULL      <- "#878787"
ARROWSIZE_ALL <- 1.2

synth_map_wk <- bind_rows(edge_summ_wk4, edge_summ_wk8, edge_summ_wk12) %>%
  group_by(from, to) %>%
  summarise(
    n_nonzero = sum(decision == "nonzero", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cat      = case_when(
      n_nonzero == 3 ~ "consistent",
      n_nonzero >= 1 ~ "sometimes",
      TRUE           ~ "insig"
    ),
    color    = ifelse(cat == "insig", COL_NULL, COL_SIG),
    style    = ifelse(cat == "insig", "dotted", "solid"),
    penwidth = case_when(
      cat == "consistent" ~ PW_CONSISTENT,
      cat == "sometimes"  ~ PW_SOMETIME,
      TRUE                ~ PW_NULL
    ),
    arrowsize = ARROWSIZE_ALL
  )

# Build DAG using identical layout code as your original
# (reuse dot_header, rank_*, horiz_*, center_* objects)

style_df_wk <- edge_summ_wk4 %>%
  select(from, to) %>%
  distinct() %>%
  left_join(synth_map_wk, by = c("from", "to"))

tg_wk <- filter(style_df_wk, from == "Temp", to == "GPP") %>% slice(1)

all_nodes_wk    <- unique(c(style_df_wk$from, style_df_wk$to,
                            top_nodes, mid2_nodes, mid3_nodes, bot_nodes))
node_lines_wk   <- vapply(all_nodes_wk,
                          function(lbl) sprintf('  "%s" [shape=ellipse, width=2.2, height=0.8];', lbl),
                          character(1))
style_no_TG_wk  <- subset(style_df_wk, !(from == "Temp" & to == "GPP"))
edge_lines_wk   <- vapply(
  seq_len(nrow(style_no_TG_wk)),
  function(i) sprintf(
    '  "%s" -> "%s" [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s];',
    style_no_TG_wk$from[i], style_no_TG_wk$to[i],
    style_no_TG_wk$color[i], style_no_TG_wk$penwidth[i],
    style_no_TG_wk$arrowsize[i], style_no_TG_wk$style[i]
  ), character(1)
)
temp_gpp_arrow_wk <- sprintf(
  '  "Temp":s -> "GPP":nw [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s, weight=70, minlen=2, constraint=false];',
  tg_wk$color, tg_wk$penwidth, tg_wk$arrowsize, tg_wk$style
)

newdot_wk <- paste(
  dot_header,
  paste(node_lines_wk, collapse = "\n"),
  rank_top, rank_mid2_redo, rank_mid3, rank_bot,
  horiz_top, horiz_mid2, horiz_mid3, keep_N_right, center_helps,
  paste(edge_lines_wk, collapse = "\n"),
  temp_gpp_arrow_wk,
  paste(center_lines, collapse = "\n"),
  "}", sep = "\n"
)

### Plot and save Fig S9 (DAGs per week and all-week summary) ###

# Solid thick  (PW=8.0) = nonzero at 20% ROPE
# Solid thin   (PW=2.5) = nonzero at 5% or 10% ROPE only
# Grey dotted  (PW=3.0) = always ambiguous/null

PW_NONZERO   <- 8.0
PW_SENSITIVE <- 2.5
PW_NULL      <- 3.0
COL_SIG      <- "#222222"
COL_NULL     <- "#878787"
ARROWSIZE_ALL <- 1.2

# Helper: classify each edge for one week 
make_style_df <- function(edge_summ_wk, week_label) {
  
  rope_wk <- rope_sensitivity_wk_pooled %>%
    filter(Week == week_label) %>%
    select(Effect, ROPE_5_pct, ROPE_10_pct, ROPE_20_pct)
  
  edge_summ_wk %>%
    mutate(Effect = paste(from, "\u2192", to)) %>%
    left_join(rope_wk, by = "Effect") %>%
    mutate(
      edge_cat = case_when(
        ROPE_20_pct <= 5                         ~ "nonzero",
        ROPE_5_pct <= 5 | ROPE_10_pct <= 5      ~ "sensitive",
        TRUE                                     ~ "insig"
      ),
      color     = if_else(edge_cat == "insig", COL_NULL, COL_SIG),
      style     = if_else(edge_cat == "insig", "dotted", "solid"),
      penwidth  = case_when(
        edge_cat == "nonzero"   ~ PW_NONZERO,
        edge_cat == "sensitive" ~ PW_SENSITIVE,
        TRUE                    ~ PW_NULL
      ),
      arrowsize = ARROWSIZE_ALL
    )
}

# customized dag builder from DiagrammeR
build_dag <- function(style_df, week_label, filename) {
  
  tg <- dplyr::filter(style_df, from == "Temp", to == "GPP") %>% dplyr::slice(1)
  # Node definitions 
  all_nodes <- unique(c(style_df$from, style_df$to,
                        top_nodes, mid2_nodes, mid3_nodes, bot_nodes))
  node_lines <- vapply(all_nodes,
                       function(lbl) sprintf('  "%s" [shape=ellipse, width=2.2, height=0.8];', lbl),
                       character(1))
  style_no_TG <- subset(style_df, !(from == "Temp" & to == "GPP"))
  edge_lines <- vapply(
    seq_len(nrow(style_no_TG)),
    function(i) sprintf(
      '  "%s" -> "%s" [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s];',
      style_no_TG$from[i], style_no_TG$to[i],
      style_no_TG$color[i], style_no_TG$penwidth[i],
      style_no_TG$arrowsize[i], style_no_TG$style[i]
    ), character(1))
  
  temp_gpp_arrow <- sprintf(
    '  "Temp":s -> "GPP":nw [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s, weight=70, minlen=2, constraint=false];',
    tg$color, tg$penwidth, tg$arrowsize, tg$style
  )
  
  newdot <- paste(
    dot_header,
    paste(node_lines,  collapse = "\n"),
    rank_top, rank_mid2_redo, rank_mid3, rank_bot,
    horiz_top, horiz_mid2, horiz_mid3, keep_N_right, center_helps,
    paste(edge_lines,  collapse = "\n"),
    temp_gpp_arrow,
    paste(center_lines, collapse = "\n"),
    "}", sep = "\n"
  )
  
  fig         <- DiagrammeR::grViz(newdot, width = "100%")
  svg_content <- DiagrammeRsvg::export_svg(fig)
  rsvg::rsvg_pdf(charToRaw(svg_content),
                 file  = filename,
                 width = 1210, height = 590)
  
  invisible(fig)
}

style_wk4  <- make_style_df(edge_summ_wk4,  "4")
style_wk8  <- make_style_df(edge_summ_wk8,  "8")
style_wk12 <- make_style_df(edge_summ_wk12, "12")

# make week-specific DAGs - Fig S9 (a-c)
fig_dag_wk4  <- build_dag(style_wk4,  "4",  "sem_dag_week4_pooled.pdf")
fig_dag_wk8  <- build_dag(style_wk8,  "8",  "sem_dag_week8_pooled.pdf")
fig_dag_wk12 <- build_dag(style_wk12, "12", "sem_dag_week12_pooled.pdf")

# check to see what they look like
fig_dag_wk4
fig_dag_wk8
fig_dag_wk12

dag_style <- bind_rows(
  edge_summ_wk4  %>% mutate(Week = "4"),
  edge_summ_wk8  %>% mutate(Week = "8"),
  edge_summ_wk12 %>% mutate(Week = "12")
) %>%
  group_by(from, to) %>%
  summarise(
    n_nonzero = sum(decision == "nonzero",  na.rm = TRUE),
    n_ambig   = sum(decision == "ambiguous", na.rm = TRUE),
    n_null    = sum(decision == "null",      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    edge_cat  = case_when(
      n_nonzero == 3                  ~ "consistent",
      n_nonzero >= 1 & n_null == 0   ~ "sometimes",
      TRUE                            ~ "insig"
    ),
    color     = if_else(edge_cat == "insig", COL_NULL, COL_SIG),
    style     = if_else(edge_cat == "insig", "dotted", "solid"),
    penwidth  = case_when(
      edge_cat == "consistent" ~ PW_NONZERO,
      edge_cat == "sometimes"  ~ PW_SENSITIVE,
      TRUE                     ~ PW_NULL
    ),
    arrowsize = ARROWSIZE_ALL
  )

# All-week summary DAG-Fig S9 (d)
dag_style_full <- edge_summ_wk4 %>%
  select(from, to) %>%
  distinct() %>%
  left_join(dag_style, by = c("from", "to"))

fig_dag_summary <- build_dag(dag_style_full, "All weeks", "sem_dag_summary_pooled.pdf")
fig_dag_summary

# Edge coefficient reference table (manually added edge labels because too complicated in DiagrammeR)
bind_rows(
  style_wk4  %>% mutate(Week = "4"),
  style_wk8  %>% mutate(Week = "8"),
  style_wk12 %>% mutate(Week = "12")
) %>%
  filter(edge_cat != "insig") %>%
  mutate(
    label     = sprintf("%.2f [%.2f, %.2f]", Mean, L89, U89),
    line_type = case_when(
      edge_cat == "nonzero"   ~ "solid thick",
      edge_cat == "sensitive" ~ "solid thin",
      TRUE                    ~ "dotted"
    )
  ) %>%
  select(Week, from, to, line_type, label) %>%
  arrange(Week, from, to) %>%
  print(n = Inf)
