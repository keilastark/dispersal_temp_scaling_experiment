# Supplementary Figures for:
# 
# "Dispersal transiently modifies the temperature dependence of ecosystem productivity after an extreme thermal fluctuation"
#
# K. Stark, Dec 2025
#
# Notes: 
# Code to reproduce Fig S4 (Directed Acyclic Graph showing SEM results) can be found in "SEM_analyses_and_DAG.R"
# Code to reproduce Fig S12a (regressions showing LMM of dispersal effects on phytoplankton size spectra) can be found in "size_spectra_analyses.R"

# Nitrate/phosphate concentrations by wattage treatment over time (Fig. S1)--------- 
nuts<-read.csv("nutrients.csv")
nuts <- nuts %>%
  mutate(dispersal = factor(dispersal),
         sample_date =)

nuts$sample_date <- as.character(nuts$sample_date)
week_mapping <- c("2018-06-27" = 2, "2018-07-11" = 4, "2018-07-25" = 6, 
                  "2018-08-08" = 8, "2018-08-22" = 10, "2018-09-05" = 12)
nuts <- nuts %>%
  mutate(week_number = dplyr::recode(date, !!!week_mapping)) 

nuts$date<-as.Date(nuts$date)

sampling_dates <- ymd(sampling_dates)       

nuts_filtered <- nuts %>% 
  rename(sample_date = date) %>%            
  filter(sample_date %in% sampling_dates) 

nuts$wattage<-  factor(nuts$wattage,levels=c(0,100,200,300))

nitrates_all_dates <- ggplot(nuts, aes(x=date,y=mean_nit,group=interaction(date,wattage))) + geom_boxplot(aes(fill=wattage))+ theme_classic()+scale_fill_brewer(palette="RdBu",direction = -1,name="Wattage") +ylim(0,7)+ylab("[Nitrate]") + xlab("") +
  theme( legend.position = "right",
         legend.justification = "top",
         legend.text     = element_text(size = 11),
         legend.title    = element_text(size = 12),
         panel.border    = element_rect(colour = "black", fill = NA),
         axis.text.x     = element_text(size = 12),
         axis.text.y     = element_text(size = 12),
         axis.title.y    = element_text(size = 13),
         axis.title.x    = element_text(size = 13),
         plot.title.position = "plot",
         plot.margin     = margin(12,12,12,12)
  )

phosphates_all_dates <- ggplot(nuts, aes(x=date,y=mean_pho,group=interaction(date,wattage)))+geom_boxplot(aes(fill=factor(wattage)))+theme_classic()+scale_fill_brewer(palette="RdBu",direction = -1,name="Wattage")  +ylim(0,2.1)+ylab("[Phosphate]")+xlab("Date")+
  theme_classic(base_size = 13) + 
  theme(legend.position = "none",
        legend.text     = element_text(size = 11),
        legend.title    = element_text(size = 12),
        panel.border    = element_rect(colour = "black", fill = NA),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        axis.title.y    = element_text(size = 13),
        axis.title.x    = element_text(size = 13),
        plot.title.position = "plot",
        plot.margin     = margin(12,12,12,12))

nutrients_combined <- nitrates_all_dates /  phosphates_all_dates + 
  plot_layout(widths = c(1,1)) +
  plot_annotation(
    
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
  ) &
  theme(
    plot.tag.position = c(-0.02, 1.07), # small positive x, consistent across panels
    plot.tag = element_text(size = 16, face = "bold", family = "Helvetica",
                            hjust = 0, vjust = 1)
  )

nutrients_combined
setwd("~/Documents/Pond_MS/figures")
ggsave(nutrients_combined, file="nutrients_combined.pdf", width=6.9,height=5.8)


#### Testing all-week version of H1 LMM (Fig. S4) ####

draws_u <- posterior::as_draws_df(fit_unified_h1)
get_u   <- function(v) {
  if (v %in% names(draws_u)) as.numeric(draws_u[[v]]) else rep(0, nrow(draws_u))
}

make_row_u <- function(draws, week, param) {
  tibble(
    Week      = factor(week, levels = rev(WLVL)),
    Parameter = param,
    Mean      = mean(draws),
    L89       = quantile(draws, 0.055),
    U89       = quantile(draws, 0.945),
    HDIinROPE = hdi_in_rope_pct(draws, thr),
    Decision  = case_when(
      HDIinROPE >= 95 ~ "null",
      HDIinROPE <=  5 ~ "nonzero",
      TRUE            ~ "ambiguous"
    )
  )
}

b_temp_u <- get_u("b_inv_kTc")

unified_rows <- purrr::map_dfr(WLVL, function(wk) {
  wk_num <- gsub("Week ", "", wk)
  
  none_d  <- b_temp_u + get_u(paste0("b_inv_kTc:week_numberWeek", wk_num, ":dispersalnone"))
  low_d   <- b_temp_u + get_u(paste0("b_inv_kTc:week_numberWeek", wk_num, ":dispersallow"))
  high_d  <- b_temp_u + get_u(paste0("b_inv_kTc:week_numberWeek", wk_num, ":dispersalhigh"))
  
  delta_low_d  <- get_u(paste0("b_inv_kTc:week_numberWeek", wk_num, ":dispersallow")) -
    get_u(paste0("b_inv_kTc:week_numberWeek", wk_num, ":dispersalnone"))
  delta_high_d <- get_u(paste0("b_inv_kTc:week_numberWeek", wk_num, ":dispersalhigh")) -
    get_u(paste0("b_inv_kTc:week_numberWeek", wk_num, ":dispersalnone"))
  
  bind_rows(
    make_row_u(none_d,       wk, "E_GPP (baseline, no dispersal)"),
    make_row_u(low_d,        wk, "E_GPP: low dispersal"),
    make_row_u(high_d,       wk, "E_GPP: high dispersal"),
    make_row_u(delta_low_d,  wk, "Temp \u00d7 Disp (low \u2212 none)"),
    make_row_u(delta_high_d, wk, "Temp \u00d7 Disp (high \u2212 none)")
  )
})

# Now panel_order_u and panels_keep_u use exactly these labels
panel_order_u <- c(
  "E_GPP (baseline, no dispersal)",
  "E_GPP: low dispersal",
  "E_GPP: high dispersal",
  "Temp \u00d7 Disp (low \u2212 none)",
  "Temp \u00d7 Disp (high \u2212 none)"
)

unified_rows <- unified_rows %>%
  mutate(
    Parameter  = factor(Parameter, levels = panel_order_u),
    color_flag = if_else(Decision == "nonzero", "nonzero", "other")
  )

panel_rename_pw <- c(
  "Temperature (1/kTc)"          = "E_GPP (baseline, no dispersal)",
  "Temp \u00d7 Dispersal (low)"  = "Temp \u00d7 Disp (low \u2212 none)",
  "Temp \u00d7 Dispersal (high)" = "Temp \u00d7 Disp (high \u2212 none)"
)

panels_keep <- c(
  "E_GPP (baseline, no dispersal)",
  "Temp \u00d7 Disp (low \u2212 none)",
  "Temp \u00d7 Disp (high \u2212 none)"
)

perweek_sub <- perweek_rows %>%
  filter(Parameter %in% c(
    "Temperature (1/kTc)",
    "Temp \u00d7 Dispersal (low)",
    "Temp \u00d7 Dispersal (high)"
  )) %>%
  mutate(Parameter = dplyr::recode(as.character(Parameter), !!!panel_rename_pw)) %>%
  mutate(Parameter = factor(Parameter, levels = panels_keep))

rope_rects_pw_sub <- tibble(
  Parameter = factor(panels_keep, levels = panels_keep),
  xmin = -thr, xmax = thr
)

fig_forest_perweek_sub <- ggplot(
  perweek_sub,
  aes(y = Week, x = Mean, colour = color_flag)
) +
  geom_rect(
    data = rope_rects_pw_sub,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "grey85", alpha = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black",
             linewidth = 0.5) +
  geom_errorbarh(aes(xmin = L89, xmax = U89), height = 0.25, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_colour_manual(
    values = c("nonzero" = col_nonzero, "other" = col_other)
  ) +
  facet_wrap(~Parameter, nrow = 1, scales = "free_x") +
  labs(y = "Week",
       x = "Posterior mean with 89% HDI (shaded band = ROPE)") +
  forest_theme



fig_s4_combined <- (fig3 + ggtitle("(a)")) /
  (fig_forest_perweek_sub + ggtitle("(b)")) +
  plot_layout(heights = c(1, 1))

fig_s4_combined




#### 
#### Prior sensitivity analysis of H1 LMM (Fig. S5) ####

# Refit per-week models under three prior specifications and compare posteriors for the key parameters: E_GPP (inv_kTc) and dispersal:temp interactions

# Define different prior sets:
# Original (informative, as in manuscript)
priors_original <- c(
  prior(normal(1.5, 5),     class = "Intercept"),
  prior(normal(0.65, 0.28), class = "b", coef = "inv_kTc"),
  prior(normal(0, 0.5),     class = "b", coef = "inv_kTc:dispersallow"),
  prior(normal(0, 0.5),     class = "b", coef = "inv_kTc:dispersalhigh"),
  prior(exponential(1),     class = "sigma"),
  prior(exponential(1),     class = "sd")
)

# Wider E_GPP prior (double SD)
priors_wide <- c(
  prior(normal(1.5, 5),     class = "Intercept"),
  prior(normal(0.65, 0.6), class = "b", coef = "inv_kTc"),
  prior(normal(0, 1),       class = "b", coef = "inv_kTc:dispersallow"),
  prior(normal(0, 1),       class = "b", coef = "inv_kTc:dispersalhigh"),
  prior(exponential(1),     class = "sigma"),
  prior(exponential(1),     class = "sd")
)

# Weakly informative (flat on E_GPP, wide on interactions)
priors_weak <- c(
  prior(normal(1.5, 5),     class = "Intercept"),
  prior(normal(0, 10),       class = "b", coef = "inv_kTc"),
  prior(normal(0, 10),       class = "b", coef = "inv_kTc:dispersallow"),
  prior(normal(0, 10),       class = "b", coef = "inv_kTc:dispersalhigh"),
  prior(exponential(1),     class = "sigma"),
  prior(exponential(1),     class = "sd")
)

prior_sets <- list(
  "Informative (manuscript)" = priors_original,
  "Wide"                     = priors_wide,
  "Uninformative"       = priors_weak
)

#re-fit for each prior category
weeks_sensitivity <- c("Week 2", "Week 4", "Week 6", "Week 8", "Week 10", "Week 12")

fit_sensitivity <- function(wk, prior_set) {
  dat_d <- dat2 %>%
    filter(week_number == wk) %>%
    tidyr::drop_na(ln_GPP_umol_gC, inv_kTc, dispersal, metacommunity)
  
  brm(
    ln_GPP_umol_gC ~ inv_kTc + inv_kTc:dispersal + (1 | metacommunity),
    data      = dat_d,
    family    = gaussian(),
    prior     = prior_set,
    chains    = 4, cores = 4, iter = 5000, warmup = 1000,
    backend   = "cmdstanr",
    control   = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    seed      = 123
  )
}

sensitivity_fits <- purrr::map(prior_sets, function(ps) {
  setNames(
    purrr::map(weeks_sensitivity, ~fit_sensitivity(.x, ps)),
    weeks_sensitivity
  )
})


extract_sensitivity <- function(fit, week, prior_name) {
  d    <- posterior::as_draws_df(fit)
  get_d <- function(v) if (v %in% names(d)) as.numeric(d[[v]]) else rep(0, nrow(d))
  
  bT   <- get_d("b_inv_kTc")
  bTL  <- get_d("b_inv_kTc:dispersallow")
  bTH  <- get_d("b_inv_kTc:dispersalhigh")
  
  summarise_par <- function(draws, param) {
    tibble(
      Prior     = prior_name,
      Week      = week,
      Parameter = param,
      Mean      = mean(draws),
      L89       = quantile(draws, 0.055),
      U89       = quantile(draws, 0.945),
      L95       = quantile(draws, 0.025),
      U95       = quantile(draws, 0.975)
    )
  }
  
  bind_rows(
    summarise_par(bT,        "E_GPP (none dispersal)"),
    summarise_par(bT + bTL,  "E_GPP (low dispersal)"),
    summarise_par(bT + bTH,  "E_GPP (high dispersal)"),
    summarise_par(bTL,       "Delta E_GPP (low - none)"),
    summarise_par(bTH,       "Delta E_GPP (high - none)")
  )
}

sensitivity_results <- purrr::map_dfr(names(sensitivity_fits), function(pname) {
  purrr::map_dfr(names(sensitivity_fits[[pname]]), function(wk) {
    extract_sensitivity(sensitivity_fits[[pname]][[wk]], wk, pname)
  })
}) %>%
  mutate(
    Prior = factor(Prior, levels = names(prior_sets)),
    Week  = factor(Week,  levels = WLVL)
  )

#make fig S5
params_show <- c(
  "E_GPP (none dispersal)",
  "E_GPP (low dispersal)",
  "E_GPP (high dispersal)",
  "Delta E_GPP (low - none)",
  "Delta E_GPP (high - none)"
)

sensitivity_plot_data <- sensitivity_results %>%
  filter(Parameter %in% params_show) %>%
  mutate(Parameter = factor(Parameter, levels = params_show))

prior_colours <- c(
  "Informative (manuscript)" = "#2c7bb6",
  "Wide"                     = "#d7191c",
  "Uninformative"       = "#1a9641"
)

prior_shapes <- c(
  "Informative (manuscript)" = 19,
  "Wide"                     = 17,
  "Uninformative"       = 15
)

fig_sensitivity <- ggplot(
  sensitivity_plot_data,
  aes(y = Prior, x = Mean, colour = Prior, shape = Prior)
) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = L89, xmax = U89),
    height = 0.25, linewidth = 0.7
  ) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = prior_colours) +
  scale_shape_manual(values = prior_shapes) +
  facet_grid(
    rows = vars(Week),
    cols = vars(Parameter),
    scales = "free_x"
  ) +
  labs(
    x      = "Posterior mean with 89% HDI",
    y      = NULL,
    colour = "Prior specification",
    shape  = "Prior specification"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text       = element_text(size = 10, family = "Helvetica"),
    panel.border     = element_rect(colour = "black", fill = NA),
    axis.text.y      = element_text(size = 9, family = "Helvetica"),
    axis.text.x      = element_text(size = 9, family = "Helvetica"),
    axis.title.x     = element_text(size = 11, family = "Helvetica",
                                    margin = margin(t = 8)),
    legend.text      = element_text(size = 10, family = "Helvetica"),
    legend.title     = element_text(size = 11, family = "Helvetica"),
    panel.spacing.x  = unit(0.4, "cm"),
    panel.spacing.y  = unit(0.3, "cm")
  )

fig_sensitivity

ggsave(fig_sensitivity,
       file  = "fig_prior_sensitivity.pdf",
       width = 11, height = 9.7)

# table version
sensitivity_table <- sensitivity_results %>%
  filter(Parameter %in% params_show) %>%
  mutate(
    across(c(Mean, L89, U89), ~round(.x, 3))
  ) %>%
  select(Week, Parameter, Prior, Mean, L89, U89) %>%
  arrange(Week, Parameter, Prior)

print(sensitivity_table, n = Inf)
##### Frequency plot of size bins (Fig. S11) ####
phyto_summary_disp <- phyto_size_long %>%
  group_by(week_number, dispersal, size_bin) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE),
            se_prop   = sd(prop, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

phyto_size_long_collapsed <- phyto_size_long %>% # combined the 3 largest (most sparse) size bins
  mutate(
    size_bin = fct_recode(
      size_bin,
      "≥ 25 µm" = "25_30 µm",
      "≥ 25 µm" = "30_35 µm",
      "≥ 25 µm" = "35_40 µm"
    )
  )

phyto_summary_disp <- phyto_size_long_collapsed %>%
  group_by(week_number, dispersal, size_bin) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    se_prop   = sd(prop,  na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

size_levels <- c("0_5 µm","5_10 µm","10_15 µm","15_20 µm","20_25 µm","≥ 25 µm")

phyto_summary_disp <- phyto_summary_disp %>%
  mutate(size_bin = factor(size_bin, levels = size_levels))

size_labels <- c(
  "0_5 µm"   = "0 ≤ 5",
  "5_10 µm"  = "5 ≤ 10",
  "10_15 µm" = "10 ≤ 15",
  "15_20 µm" = "15 ≤ 20",
  "20_25 µm" = "20 ≤ 25",
  "≥ 25 µm"  = "≥ 25"
)

disp_colors <- c("none"="#003f5c","low"="#ff6361","high"="#ffa600")

disp_freq_size_plot <- ggplot(
  phyto_summary_disp,
  aes(x = size_bin, y = mean_prop, fill = dispersal)
) +
  geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.9) +
  geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
                position = position_dodge(width = 0.8), width = 0.3) +
  facet_wrap(
    ~ week_number, nrow = 2,
    labeller = labeller(week_number = c("2"="Week 2","4"="Week 4","6"="Week 6",
                                        "8"="Week 8","10"="Week 10","12"="Week 12"))
  ) +
  scale_x_discrete(breaks = size_levels, labels = size_labels) +
  scale_fill_manual(values = disp_colors, name = "Dispersal") +
  labs(
    x = "Phytoplankton size bin (Equivalent Spherical Diameter, µm)",
    y = "Mean proportion (±SE)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.justification="top",
    
    legend.margin = margin(0, 0, 0, 2),   
    plot.margin   = margin(5, 5, 5, 5) ,
    legend.title = element_text(size = 12, family = "Helvetica"),
    strip.text    = element_text(size = 12, family = "Helvetica"),
    legend.text   = element_text(size = 11, family = "Helvetica"),
    axis.title    = element_text(size = 14, family = "Helvetica"),
    axis.title.x = element_text(size=14,vjust=-1),
    axis.text     = element_text(size = 11, family = "Helvetica"),
    axis.title.y  = element_text(vjust = 3),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    panel.border  = element_rect(color = "black", fill = NA)
  )

#ggsave(disp_freq_size_plot,file="disp_freq_size_plot.pdf",width=6.1,height=4.8) # Fig. S11


# Mean phyto body size by wattage treatment over time (Fig. S12b) ----------

vol_sum <- exp_data %>%
  mutate(week_number = factor(week_number, levels = c("2","4","6","8","10","12"))) %>%
  group_by(week_number, wattage) %>%
  summarise(
    mean_vol  = mean(volume_mean, na.rm = TRUE),
    se_vol    = sd(volume_mean,  na.rm = TRUE) / sqrt(sum(!is.na(volume_mean))),
    mean_temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  )

shape_vals <- c("0"=21, "100"=24, "200"=22, "300"=23) 
rdbu_cont <- rev(brewer.pal(11, "RdBu"))

phyto_size_watt_time <- ggplot(vol_sum, aes(x = week_number, y = mean_vol,  fill = mean_temp,shape = wattage)) +
  geom_errorbar(aes(ymin = mean_vol - se_vol, ymax = mean_vol + se_vol),
                position = position_dodge(width = 0.5) , width = 0.5, linewidth = 0.4, color = "black") +
  geom_point(position = position_dodge(width = 0.5) , size = 3, color = "black") +  
  scale_fill_gradientn(colors = rdbu_cont, name = "Mean temp (°C)",
                       breaks = c(20, 24, 28)) +
  scale_shape_manual(values = shape_vals, name = "Wattage") +
  labs(
    x = "Week",
    y = expression("Mean phytoplankton cell volume (µm"^3*")")
  ) +
  theme_classic(base_size = 13) + guides(
    fill  = guide_colorbar(order = 1, title.position = "top"),
    shape = guide_legend(order = 2,  title.position = "top")
  )+
  theme(
    legend.position = "right",     
    legend.justification = "top",
    legend.box = "vertical",           
    legend.box.just = "left",
    legend.title = element_text(size = 12, family = "Helvetica"),
    strip.text  = element_text(size = 12, family = "Helvetica"),
    legend.text = element_text(size = 11, family = "Helvetica"),
    axis.title  = element_text(size = 14, family = "Helvetica"),
    axis.text   = element_text(size = 12, family = "Helvetica"),
    axis.title.x = element_text(vjust = -0.3),
    axis.title.y = element_text(vjust =  3),
    panel.border = element_rect(color = "black", fill = NA)
  )  

# Note: need to have run size_spectra_analyses.R to get p_week_models ggplot object
combined_size <- p_week_models  + plot_spacer()+ phyto_size_watt_time + 
  plot_layout(widths = c(1.7,0.1,0.9)) +
  plot_annotation(
    
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(plot.margin = margin(t = 25, r = 6, b = 6, l = 16))  
  ) &
  theme(
    plot.tag.position = c(-0.02, 1.07), 
    plot.tag = element_text(size = 16, face = "bold", family = "Helvetica",
                            hjust = 0, vjust = 1)
  )
combined_size 

# Zooplankton abundance vs inedible phytolankton abundance (Fig. S12c)

phyto_size_v_zoop <-  exp_data %>%
  filter(prop_inedible >0) %>%
  ggplot( aes(x = log(total_zoop_abundance), y = log(prop_inedible), color = temp)) +
  geom_point(size = 3, stroke=0,alpha=0.8) +  # no outlines
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_gradientn(colors = rdbu_cont, name = "Temp (°C)",
                        breaks = c(20, 24, 28)) +
  labs(
    x = expression("ln Zooplankton biomass (µg L"^-1*")"),
    y = expression("ln large phyto density (ESD >30µm, ind L"^-1*")")
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",                      # suppress legends
    legend.title = element_text(size = 13, family = "Helvetica"),
    strip.text    = element_text(size = 12, family = "Helvetica"),
    legend.text   = element_text(size = 10, family = "Helvetica"),
    axis.title    = element_text(size = 14, family = "Helvetica"),
    axis.text     = element_text(size = 12, family = "Helvetica"),
    axis.title.x  = element_text(vjust = -2),
    axis.title.y  = element_text(vjust =  3),
    panel.border  = element_rect(color = "black", fill = NA)
  )+
  theme(
    legend.key.size = unit(0.4, "cm"),          
    legend.text = element_text(size = 9),       
    legend.title = element_text(size = 10),
    legend.spacing.y = unit(0.1, "cm"),         
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),plot.margin = margin(t = 6, r = 6, b = 6, l = 16)
  ) 

#ggsave(phyto_size_v_zoop, file="phyto_zoop_size.pdf",width=4,height=3.4)


### Testing for week-specific differences in relative abundance of Oocystis sp. (Fig S15) ###
library(dplyr)
library(tidyr)
library(lme4)
library(brms)
library(ggplot2)

#  Merge phyto count data with experimental metadata as in the main phyto diversity script

Phyto_data<-read.csv("phytoplankton_diversity_data.csv")
Phyto_data <- Phyto_data %>%
  mutate(tank = as.factor(tank))

meta_link <- exp_data %>%
  dplyr::select(tank, week_number, dispersal, metacommunity, wattage) %>%
  distinct()

date_to_week <- tibble(
  sample_date = sort(unique(Phyto_data$sample_date)),
  week_number = c(4L, 8L, 12L)  
)

print(date_to_week)
# should be
# 2018-07-11 -> 4
# 2018-08-08 -> 8
# 2018-09-05 -> 12

Phyto_labelled <- Phyto_data %>%
  left_join(date_to_week, by = "sample_date") %>%
  left_join(meta_link, by = c("tank", "week_number"))

# calculate Oocystis relative abundance

Oocystis_ra <- Phyto_labelled %>%
  mutate(
    # Combine large and small Oocystis as one genus
    oocystis_count = rowSums(
      dplyr::select(., Large.oocystis, Small.oocystis),
      na.rm = TRUE
    ),
    total_count = rowSums(
      dplyr::select(., Ceratium:Cryptomonas.sp.),
      na.rm = TRUE
    ),
    oocystis_prop = oocystis_count / total_count,
    # Logit transform for modelling (add small constant to handle 0s)
    oocystis_logit = car::logit(oocystis_prop, adjust = 0.025)
  ) %>%
  filter(!is.na(dispersal), total_count > 0) %>%
  mutate(
    dispersal     = factor(dispersal, levels = c("none", "low", "high")),
    metacommunity = factor(metacommunity),
    week_number   = factor(week_number),
    tank          = factor(tank)
  )

Oocystis_ra %>%
  group_by(dispersal, week_number) %>%
  summarise(
    mean_prop = mean(oocystis_prop, na.rm = TRUE),
    sd_prop   = sd(oocystis_prop,   na.rm = TRUE),
    n         = n(),
    .groups   = "drop"
  ) %>%
  print(n = Inf)

# linear mixed model on logit-transformed relative abundance (See Supp Methods)

lmm_oocystis <- lme4::lmer(
  oocystis_logit ~ dispersal * week_number + wattage +
    (1 | metacommunity),
  data = Oocystis_ra,
  REML = TRUE
)

summary(lmm_oocystis)

# Likelihood ratio test for dispersal effect
lmm_null <- lme4::lmer(
  oocystis_logit ~ week_number + wattage + (1 | metacommunity),
  data = Oocystis_ra, REML = FALSE
)
lmm_full <- lme4::lmer(
  oocystis_logit ~ dispersal * week_number + wattage + (1 | metacommunity),
  data = Oocystis_ra, REML = FALSE
)
anova(lmm_null, lmm_full)

# Bayesian version (consistent with rest of manuscript)

priors_oocystis <- c(
  prior(normal(0, 1),   class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(1), class = "sigma"),
  prior(exponential(1), class = "sd")
)

fit_oocystis <- brm(
  oocystis_logit ~ dispersal * week_number + wattage +
    (1 | metacommunity),
  data    = Oocystis_ra,
  family  = gaussian(),
  prior   = priors_oocystis,
  chains  = 4, cores = 4, iter = 4000, warmup = 1000,
  backend = "cmdstanr",
  seed    = 123
)

summary(fit_oocystis)

# extract posterior dispersal effect within each week
draws_ooc <- posterior::as_draws_df(fit_oocystis)
get_ooc   <- function(v) if (v %in% names(draws_ooc)) as.numeric(draws_ooc[[v]]) else rep(0, nrow(draws_ooc))

b_int       <- get_ooc("b_Intercept")
b_low       <- get_ooc("b_dispersallow")
b_high      <- get_ooc("b_dispersalhigh")
b_wk8       <- get_ooc("b_week_number8")
b_wk12      <- get_ooc("b_week_number12")
b_low_wk8   <- get_ooc("b_dispersallow:week_number8")
b_high_wk8  <- get_ooc("b_dispersalhigh:week_number8")
b_low_wk12  <- get_ooc("b_dispersallow:week_number12")
b_high_wk12 <- get_ooc("b_dispersalhigh:week_number12")

sd_ooc <- sd(Oocystis_ra$oocystis_logit, na.rm = TRUE)
thr_ooc <- 0.20 * sd_ooc

summarise_ooc <- function(draws, label) {
  tibble(
    Contrast  = label,
    Mean      = mean(draws),
    L89       = quantile(draws, 0.055),
    U89       = quantile(draws, 0.945),
    HDIinROPE = hdi_in_rope_pct(draws, thr_ooc),
    Decision  = case_when(
      HDIinROPE >= 95 ~ "null",
      HDIinROPE <=  5 ~ "nonzero",
      TRUE            ~ "ambiguous"
    )
  )
}

oocystis_contrasts <- bind_rows(
  # Week 4 (reference week)
  summarise_ooc(b_low,  "Low vs none, Week 4"),
  summarise_ooc(b_high, "High vs none, Week 4"),
  # Week 8
  summarise_ooc(b_low  + b_low_wk8,  "Low vs none, Week 8"),
  summarise_ooc(b_high + b_high_wk8, "High vs none, Week 8"),
  # Week 12
  summarise_ooc(b_low  + b_low_wk12,  "Low vs none, Week 12"),
  summarise_ooc(b_high + b_high_wk12, "High vs none, Week 12")
)

print(oocystis_contrasts, n = Inf)

## Plot Oocystis relative abundance by dispersal and week (Fig S15)

fig_oocystis <- ggplot(
  Oocystis_ra,
  aes(x = dispersal, y = oocystis_prop,
      colour = dispersal, fill = dispersal)
) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  facet_wrap(~week_number, nrow = 1,
             labeller = labeller(week_number = c(
               "4" = "Week 4", "8" = "Week 8", "12" = "Week 12"
             ))) +
  scale_colour_manual(
    values = c("none" = "#003f5c", "low" = "#ff6361", "high" = "#ffa600"),
    labels = c("None", "Low", "High"), guide = "none"
  ) +
  scale_fill_manual(
    values = c("none" = "#003f5c", "low" = "#ff6361", "high" = "#ffa600"),
    labels = c("None", "Low", "High"), guide = "none"
  ) +
  scale_x_discrete(labels = c("none" = "None", "low" = "Low", "high" = "High")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Dispersal treatment",
    y = "Oocystis relative abundance"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text       = element_text(size = 13, family = "Helvetica"),
    panel.border     = element_rect(colour = "black", fill = NA),
    axis.title       = element_text(size = 13, family = "Helvetica"),
    axis.text        = element_text(size = 11, family = "Helvetica")
  )

fig_oocystis

ggsave(fig_oocystis,
       file  = "fig_oocystis_relabun.pdf",
       width = 8, height = 3)

oocystis_plot_dat <- oocystis_contrasts %>%
  mutate(
    # Split contrast label into week and dispersal for faceting
    Week      = stringr::str_extract(Contrast, "Week \\d+"),
    Dispersal = if_else(grepl("High", Contrast), "High", "Low"),
    Dispersal = factor(Dispersal, levels = c("Low", "High")),
    Week      = factor(Week, levels = c("Week 4", "Week 8", "Week 12")),
    color_flag = if_else(Decision == "nonzero", "nonzero", "other")
  )

sd_ooc <- sd(Oocystis_ra$oocystis_logit, na.rm = TRUE)
thr_ooc <- 0.20 * sd_ooc

# Fig 16B (showing posterior HDI within ROPE = 20%)
rope_rect_ooc <- data.frame(xmin = -thr_ooc, xmax = thr_ooc)

fig_oocystis_forest <- ggplot(
  oocystis_plot_dat,
  aes(y = Dispersal, x = Mean, colour = color_flag)
) +
  geom_rect(
    data        = rope_rect_ooc,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey85", alpha = 0.7
  ) +
  geom_vline(
    xintercept = 0, linetype = "dashed",
    colour = "black", linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = L89, xmax = U89),
    height = 0.2, linewidth = 0.7
  ) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c("nonzero" = "#C62828", "other" = "grey50"),
    guide  = "none"
  ) +
  scale_y_discrete(limits = rev) +
  facet_wrap(~Week, nrow = 1) +
  labs(
    x = "Dispersal effect on relative abundance vs. no dispersal",
    y = ""
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text       = element_text(size = 13, family = "Helvetica"),
    panel.border     = element_rect(colour = "black", fill = NA),
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.title       = element_text(size = 13, family = "Helvetica"),
    axis.text        = element_text(size = 11, family = "Helvetica"),
    plot.caption     = element_text(size = 9,  family = "Helvetica",
                                    colour = "grey40")
  )

fig_oocystis_forest

# Stack 15a and b as in Supp Mat
fig_oocystis_combined <- fig_oocystis / fig_oocystis_forest +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = "a")

#ggsave(fig_oocystis_combined,
 #      file  = "fig_oocystis_combined.pdf",
      # width = 7.2, height = 5.5)

#### Stacked barplot of phytoplankton taxa by week and treatment (Fig S16)
phyto_clean <- Phyto_data %>%
  mutate(sample_date = as.Date(sample_date)) %>%
  select(-`Rhodomonas.sp.`, -Unknown_chlorophyte) %>%
  mutate(
    Oocystis = coalesce(Oocystis.large, 0L) + coalesce(Oocystis.small, 0L)
  ) %>%
  select(-Oocystis.large, -Oocystis.small)

exp_keys <- exp_data %>%
  mutate(sample_date = as.Date(sample_date)) %>%
  select(tank, sample_date, dispersal, week_number) %>%
  distinct()

phyto_long <- phyto_clean %>%
  left_join(exp_keys, by = c("tank", "sample_date")) %>%
  pivot_longer(
    cols = -c(sample_date, tank, dispersal, week_number),
    names_to = "species",
    values_to = "count"
  ) %>%
  mutate(count = coalesce(count, 0L))

plot_df <- phyto_long %>%
  group_by(week_number, dispersal, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(week_number, dispersal) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  mutate(
    dispersal = factor(dispersal, levels = c("none", "low", "high")),
    week_number = as.integer(week_number)
  )

phyto_rel_abun<- ggplot(plot_df, aes(x = dispersal, y = rel_abund, fill = species)) +
  geom_col(width = 0.85) +
  facet_wrap(~ factor(week_number), nrow = 1) +
  labs(x = "Dispersal by week number", y = "Relative abundance", fill = "Taxon") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +scale_fill_manual(values=c(
    "#D86020FF",
    "#F89040FF",
    "#903000FF",
    "#184068FF",
    "#F8D068FF",
    "#207890FF",
    "#C03020FF",
    "#F8C060FF",
    "#F8A058FF",
    "#F8E098FF",
    "#E85040FF",
    "#F8F8F8FF",
    "#D0D0D0FF",
    "#28A8B8FF",
    "#5C90D2FF"
  )
  )+
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 13, family = "Helvetica"),
    axis.title.x     = element_text(size = 13, family = "Helvetica",
                                    margin = margin(t = 11)),
    axis.title.y     = element_text(size = 13, family = "Helvetica",
                                    margin = margin(r = 11)),
    axis.text.x      = element_text(size = 10, family = "Helvetica"),
    axis.text.y      = element_text(size = 10, family = "Helvetica"),
    legend.text      = element_text(size = 10, family = "Helvetica"),
    legend.title     = element_text(size = 11, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(0.8, "cm")
  )


