library(vegan)
library(tidyverse)
library(multcompView)
library(patchwork)
library(ggpubr)
library(ggh4x)
library(paletteer)
library(gridExtra)
library(cowplot)
library(car)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

exp_data<-read.csv("metacom_exp_data_final.csv")
# make sure factors properly set
exp_data$metacommunity <- factor(exp_data$metacommunity)
exp_data$wattage <- factor(exp_data$wattage, levels=c("0","100","200","300"))
exp_data$dispersal <- factor(exp_data$dispersal, levels=c("none","low","high"))
exp_data$tank <- factor(exp_data$tank)
exp_data$sample_date <- as.character(exp_data$sample_date)

# Before the analyses testing hypothesis 2, here is code reproducing main manuscript Figure 2: time series of raw ln GPP and biomass by wattage and dispersal treatments

dispersal_labels <- c(
  none = "No dispersal",
  low  = "Low dispersal",
  high = "High dispersal"
)

gpp_v_time <- exp_data %>%
  ggplot(aes(x = week_number, y = ln_GPP_umol_L_h, color = wattage, group = wattage, shape = wattage)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 0.6) +
  facet_wrap(~dispersal, labeller = labeller(dispersal = dispersal_labels)) +
  scale_color_brewer(name = "Wattage treatment", palette = "RdBu", direction = -1) +
  scale_shape_manual(name = "Wattage treatment",values = c(16, 17, 15, 18)) + 
  theme_classic() +
  labs(
    x = "",
    y = expression(ln~GPP~"(µmol L h"^"-1"*")")
  ) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    strip.text = element_blank(),
    axis.title.x     = element_text(size = 14, family = "Helvetica",
                                    margin = margin(t = 11)),
    axis.title.y     = element_text(size = 14, family = "Helvetica",
                                    margin = margin(r = 11)),
    axis.text.x      = element_text(size = 13, family = "Helvetica"),
    axis.text.y      = element_text(size = 12, family = "Helvetica"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 13, family = "Helvetica"),
    legend.title     = element_text(size = 14, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(1, "cm")
  )
gpp_v_time

biom_v_time <- exp_data %>%
  ggplot(aes(x = week_number, y = ln_fluoro_ug_C_L, color = wattage, group = wattage, shape = wattage)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 0.6) +
  facet_wrap(~dispersal, labeller = labeller(dispersal = dispersal_labels)) +
  scale_color_brewer(name = "Wattage treatment", palette = "RdBu", direction = -1) +
  scale_shape_manual(name = "Wattage treatment",values = c(16, 17, 15, 18)) + 
  theme_classic() + ylim(2.3,3.3)+
  labs(
    x = "Week",
    y = expression(ln~Autotroph~biomass~"(µg C L"^"-1"*")")
  ) + scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::number_format(accuracy = 0.1)  
  ) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    strip.text = element_blank(),
    axis.title.x     = element_text(size = 14, family = "Helvetica",
                                    margin = margin(t = 11)),
    axis.title.y     = element_text(size = 14, family = "Helvetica",
                                    margin = margin(r = 11),,hjust=1),
    axis.text.x      = element_text(size = 13, family = "Helvetica"),
    axis.text.y      = element_text(size = 12, family = "Helvetica"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 13, family = "Helvetica"),
    legend.title     = element_text(size = 14, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(1, "cm")
  )

biomass_nolegend <- biom_v_time + theme(legend.position = "none",
                                        plot.margin  = margin(t = 0, r = 5, b = 2, l = 5))

gpp_nolegend <- gpp_v_time +
  labs(x = "Week") +  
  theme(
    legend.position = "none",
    axis.text.x  = element_text(color = NA),   
    axis.title.x = element_text(color = NA,    
                                margin = margin(t = 11)),
    plot.margin  = margin(t = 5, r = 5, b = 0, l = 5)
  )
# make legend
legend_gpp <- ggpubr::get_legend(
  gpp_v_time + theme(legend.position = "bottom")
)

fig2main <- plot_grid(
  gpp_nolegend,biomass_nolegend,
  labels=c("(a)","(b)"),label_size = 15, vjust=-0.5,hjust=1.3,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1)
)

# add legend to stacked fig
fig2 <- plot_grid(
  fig2main,legend_gpp,labels=NULL,
  ncol = 1,
  rel_heights = c(1,0.1)
)+
  theme(plot.margin = margin(t = 25, r = 0, b = 0, l = 30))

fig2

#ggsave(fig2,file="biom_gpp_timeseries.pdf",width=7.5,height=6.6)

###### Levene's test of residual variance by dispersal treatment ###
# Seeing if dispersal is associated with differences in variance in GPP and biomass 
# uses car package

w2 <- exp_data %>% filter(week_number == "2")
w4 <- exp_data %>% filter(week_number == "4")
w6 <- exp_data %>% filter(week_number == "6")
w8 <- exp_data %>% filter(week_number == "8")
w8n <- exp_data %>% filter(week_number == "8",dispersal == "none")
w8l <- exp_data %>% filter(week_number == "8",dispersal == "low")
w8h <- exp_data %>% filter(week_number == "8",dispersal == "high")
w10 <- exp_data %>% filter(week_number == "10")
w12 <- exp_data %>% filter(week_number == "12")

fit2_bio <- lm(ln_fluoro_ug_C_L ~ dispersal, data = w2)
fit2_gpp <- lm(ln_GPP_umol_L_h  ~ dispersal, data = w2)
leveneTest(fit2_bio,  center = mean)
leveneTest(fit2_gpp,  center = mean)

fit4_bio <- lm(ln_fluoro_ug_C_L ~ dispersal, data = w4)
fit4_gpp <- lm(ln_GPP_umol_L_h  ~ dispersal, data = w4)
leveneTest(fit4_bio,   center = mean)
leveneTest(fit4_gpp,  center = mean)

fit6_bio <- lm(ln_fluoro_ug_C_L ~ dispersal, data = w6)
fit6_gpp <- lm(ln_GPP_umol_L_h  ~ dispersal, data = w6)
leveneTest(fit6_bio,   center = mean)
leveneTest(fit6_gpp,  center = mean)

fit8_bio <- lm(ln_fluoro_ug_C_L ~ dispersal, data = w8)
leveneTest(fit8_bio, center = mean)
fit8_gpp <- lm(ln_GPP_umol_L_h  ~ dispersal, data = w8)
leveneTest(fit8_gpp, center = mean)

fit10_bio <- lm(ln_fluoro_ug_C_L ~ dispersal, data = w10)
fit10_gpp <- lm(ln_GPP_umol_L_h  ~ dispersal, data = w10)
leveneTest(fit10_bio,   center = mean)
leveneTest(fit10_gpp,  center = mean)

fit12_bio <- lm(ln_fluoro_ug_C_L ~ dispersal, data = w12)
fit12_gpp <- lm(ln_GPP_umol_L_h  ~ dispersal, data = w12)
leveneTest(fit12_bio, center = mean)
leveneTest(fit12_gpp, center = mean)


#############
#### Figure 4a (within-metacommunity coeff of variation) 
### Testing Hypothesis 2a-c (spatial insurance hypothesis)for all 3 mediator variables (biomass, size spectra, thermal traits)

metacom_sum <- exp_data %>%
  group_by(metacommunity, week_number) %>%
  summarise(
    mean_biomass        = mean(ln_fluoro_ug_C_L, na.rm = TRUE),
    sd_biomass          = sd(ln_fluoro_ug_C_L, na.rm = TRUE),
    cv_biomass          = sd_biomass / mean_biomass,
    
    mean_size_spec = mean(abs(slope_size_spectra), na.rm = TRUE),
    sd_size_spec   = sd(abs(slope_size_spectra), na.rm = TRUE),
    cv_size_spec        = sd_size_spec / mean_size_spec,
    
    mean_thermal_trait  = mean(CTI_raw, na.rm = TRUE),
    sd_thermal_trait    = sd(CTI_raw, na.rm = TRUE),
    cv_thermal_trait    = sd_thermal_trait / mean_thermal_trait,
    
    .groups = "drop"
  )

slopes_for_join <- slopes_by_week_metacommunity %>%
  ungroup() %>%
  mutate(
    week_number = as.integer(str_extract(as.character(week_number), "\\d+"))
  ) %>%
  select(metacommunity, week_number, slope_mean, slope_lo, slope_hi)

metacom_sum2 <- metacom_sum %>%
  mutate(week_number = as.character(week_number)) %>%
  left_join(
    slopes_for_join %>%
      mutate(week_number = as.character(week_number)),
    by = c("metacommunity", "week_number")
  )
meta_dispersal <- exp_data %>%
  select(metacommunity, dispersal) %>%
  distinct()

metacom_sum3 <- metacom_sum2 %>%
  left_join(meta_dispersal, by = "metacommunity")

week_labs <- c(`2` = "Week 2",`4` = "Week 4",`6` = "Week 6",`8` = "Week 8",`10` = "Week 10", `12` = "Week 12")

cv_long <- metacom_sum3 %>%
  filter(week_number %in% c(2, 4, 6, 8, 10, 12)) %>%
  mutate(week_number = factor(week_number, levels = c(2, 4, 6, 8, 10, 12))) %>%
  select(
    metacommunity,
    dispersal,
    week_number,
    cv_biomass,
    cv_size_spec,
    cv_thermal_trait
  ) %>%
  pivot_longer(
    cols      = starts_with("cv_"),
    names_to  = "trait",
    values_to = "cv"
  ) %>%
  mutate(
    trait = case_when(
      trait == "cv_biomass"        ~ "Biomass CV",
      trait == "cv_size_spec"      ~ "Size spectra CV",
      trait == "cv_thermal_trait"  ~ "Thermal trait index CV",
      TRUE ~ trait
    ),
    trait = factor(trait,
                   levels = c("Biomass CV", "Size spectra CV", "Thermal trait index CV"))
  )

# Compute mean + SE 
cv_summary <- cv_long %>%
  group_by(week_number, dispersal, trait) %>%
  summarise(
    mean_cv = mean(cv, na.rm = TRUE),
    se_cv   = sd(cv, na.rm = TRUE) / sqrt(sum(!is.na(cv))),
    .groups = "drop"
  )

# plot as in Fig 4a (biomass and size spectra) and Fig 5b (thermal trait composition)
cv_summary_plot <- ggplot(
  cv_summary,
  aes(
    x     = week_number,
    y     = mean_cv,
    color = dispersal,
    group = dispersal
  )
) +geom_jitter(
  data = cv_summary,
  aes(x = week_number, y = mean_cv, color = dispersal),
  width = 0.15,
  alpha = 0.25,
  size  = 1.5
) + 
  geom_point(position = position_dodge(width = 0.5),
             size = 2.5) +
  geom_errorbar(
    aes(ymin = mean_cv - se_cv, ymax = mean_cv + se_cv),
    width = 0.12,
    position = position_dodge(width = 0.5),
    linewidth = 0.6
  ) +
  facet_wrap(~ trait, scales = "free_y") +
  scale_color_manual(values = disp_pal) +
  labs(
    x = "Week",
    y = "Within-metacommunity coefficient of variation",
    color = "Dispersal"
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border     = element_rect(color = "black", fill = NA),
    strip.text       = element_text(size = 13),
    axis.text        = element_text(size = 11),
    axis.title       = element_text(size = 13),
    legend.position  = "bottom"
  )

cv_summary_plot
run_anova_tukey <- function(df, trait, week) {
  
  df_w <- df %>% filter(week_number == week)
  
  fmla <- as.formula(paste(trait, "~ dispersal"))
  
  model <- aov(fmla, data = df_w)
  
  # ANOVA 
  anova_tab <- broom::tidy(model) %>%
    filter(term == "dispersal") %>%
    select(statistic, p.value) %>%
    rename(F_value = statistic, anova_p = p.value)
  
  # Tukey posthoc testing for which weeks CV varied by dispersal treatment
  tuk <- pairs(emmeans(model, ~ dispersal), adjust = "tukey")
  tuk_tab <- broom::tidy(tuk) %>%
    select(contrast, estimate, p.value) %>%
    rename(tukey_p = p.value)
  tuk_tab %>%
    mutate(
      trait = trait,
      week  = week,
      F_value = anova_tab$F_value,
      anova_p = anova_tab$anova_p,
      sig = ifelse(anova_p < 0.05, "*", "ns")
    ) %>%
    select(trait, week, F_value, anova_p, sig, everything())
}

response_list <- list(
  cv_biomass       = c(2,4,6,8,10,12),
  cv_size_spec     = c(2,4,6,8,10,12),
  cv_thermal_trait = c(4,8,12)
)

anova_results_all <- imap_dfr(
  response_list,
  ~ map_dfr(.y, ~ run_anova_tukey(metacom_sum3, response = .x, week = .x))
)

anova_results_clean <- anova_results_all %>%
  mutate(
    response = recode(
      response,
      cv_biomass       = "Biomass CV",
      cv_size_spec     = "Size spectra CV",
      cv_thermal_trait = "Thermal trait index CV"
    )
  ) %>%
  arrange(response, week, contrast)

# nice labels for week facet 
week_labs <- c(
  `2`  = "Week 2",
  `4`  = "Week 4",
  `6`  = "Week 6",
  `8`  = "Week 8",
  `10` = "Week 10",
  `12` = "Week 12"
)

cv_long <- metacom_sum3 %>%
  pivot_longer(
    cols      = starts_with("cv_"),
    names_to  = "trait",
    values_to = "cv"
  ) %>%
  mutate(
    trait = dplyr::recode(
      trait,
      "cv_biomass"       = "Biomass CV",
      "cv_size_spec"     = "Size spectra CV",
      "cv_thermal_trait" = "Thermal trait index CV"
    ),
    trait = factor(
      trait,
      levels = c("Biomass CV", "Size spectra CV", "Thermal trait index CV")
    ),
    week_number = factor(
      week_number,
      levels = c(2, 4, 6, 8, 10, 12)
    )
  )


cv_summary_long <- cv_long %>%
  group_by(trait, week_number, dispersal) %>%
  summarise(
    mean_cv = mean(cv, na.rm = TRUE),
    se_cv   = sd(cv, na.rm = TRUE) / sqrt(sum(!is.na(cv))),
    .groups = "drop"
  )

cv_summary_plot <- cv_summary_long %>% 
  filter(trait %in% c("Biomass CV","Size spectra CV")) %>%
  ggplot(aes(
    x        = week_number,
    y        = mean_cv,
    group    = dispersal
  )
  ) +
  geom_point(
    aes(color = dispersal, shape = dispersal),
    position  = position_dodge(width = 0.4),
    size      = 3
  ) +
  geom_errorbar(
    aes(
      ymin = mean_cv - se_cv,
      ymax = mean_cv + se_cv,
      color = dispersal
    ),
    width     = 0.5,
    alpha     = 0.8,
    position  = position_dodge(width = 0.4)
  ) + facet_wrap(~ trait, scales = "free_y") +
  scale_color_manual(values = disp_pal) +
  scale_fill_manual(values = disp_pal) +  
  labs(
    x = "Week",
    y = "Coeff. var. (within-metacommunity)",
    color = "Dispersal",
    shape = "Dispersal",
    fill = "Dispersal"
  ) + scale_y_continuous(
    breaks = function(lims) seq(lims[1], lims[2], length.out = 3),
    labels = function(x) formatC(x, digits = 1, format = "fg")
  ) + 
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 12, family = "Helvetica"),
    axis.title.x     = element_text(size = 13, family = "Helvetica",
                                    margin = margin(t = 11)),
    axis.title.y     = element_text(size = 13, family = "Helvetica",
                                    margin = margin(r = 11)),
    axis.text.x      = element_text(size = 10, family = "Helvetica"),
    axis.text.y      = element_text(size = 10, family = "Helvetica"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11, family = "Helvetica"),
    legend.title     = element_text(size = 12, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(0.8, "cm")
  )
cv_summary_plot # as in Figure 4a

#ggsave(cv_summary_plot,file="cv_summary_plot.pdf",width=5.8,height=3.4)


# HYPOTHESIS 2B: Size spectra ---------------------------------------------
library(RColorBrewer)
# wrangle size bin data 
phyto_mids <- exp_data %>%
  mutate(sample_id = row_number()) %>%
  pivot_longer(
    cols = dplyr::starts_with("X"),
    names_to  = "bin",
    values_to = "count"
  ) %>%
  # "X0_5um" → min_size = 0, max_size = 5
  mutate(
    bin      = str_remove(bin, "^X"),
    bin      = str_remove(bin, "um$"),
    bin      = str_replace(bin, "_", "-")
  ) %>%
  separate(bin, into = c("min_size", "max_size"),
           sep = "-", convert = TRUE) %>%
  mutate(
    size_mid = (min_size + max_size) / 2,  
    count    = round(count)                
  )

phyto_size_spectra <- phyto_mids %>%
  filter(count > 0) %>%  
  mutate(
    log_size  = log10(size_mid),
    log_count = log10(count),
    week_f    = factor(week_number)
  )

# by wattage
phyto_size_spectra<-phyto_size_spectra %>% mutate(
  week_f = factor(
    week_f,
    levels = c("2","4","6","8","10","12"),
    labels = c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12"),
    ordered = TRUE
  ))

# mutate to get size bin (um) info
phyto_size_spectra <- exp_data %>%
  select(
    tank, week_number, wattage, dispersal, temp,metacommunity,
    X0_5um:X35_40um  
  ) %>%
  pivot_longer(
    cols = X0_5um:X35_40um,
    names_to  = "bin",
    values_to = "count"
  ) %>%
  mutate(
    bin      = str_remove(bin, "^X"),
    bin      = str_remove(bin, "um$"),
    bin      = str_replace(bin, "_", "-")
  ) %>%
  separate(bin, into = c("min_size", "max_size"),
           sep = "-", convert = TRUE) %>%
  mutate(
    size_mid  = (min_size + max_size) / 2,
    count     = round(count),
    log_size  = log10(size_mid),
    log_count = log10(count),
    week_f    = factor(week_number)
  ) %>%
  filter(is.finite(log_count), is.finite(log_size))


### Fitting Bayesian LMM for size spectra (hypothesis 2b test) ##
#  nest data by week and fit separate brms model for each
phyto_brm_dat <- phyto_size_spectra %>%
  mutate(week_f = factor(week_number))

weeks_nested <- phyto_brm_dat %>%
  group_by(week_f) %>%
  nest()

weeks_models <- weeks_nested %>%
  mutate(
    model = map(
      data,
      ~ brm(
        log_count ~ log_size+log_size:dispersal,
        data    = .x,
        family  = gaussian(),
        cores   = 4,
        chains  = 4,
        iter    = 4000,
        control = list(adapt_delta = 0.95)
      )
    )
  )

coef_by_week <- map2_dfr(
  weeks_models$model,
  weeks_models$week_f,
  ~ {
    fe <- as.data.frame(fixef(.x))  
    fe$term <- rownames(fe)
    
    tibble(
      week_f   = .y,
      term     = fe$term,
      estimate = fe$Estimate,
      lower95  = fe$Q2.5,
      upper95  = fe$Q97.5
    )
  }
) %>%
  arrange(as.numeric(as.character(week_f)), term)

coef_by_week

sdy <- sd(phyto_brm_dat$log_count, na.rm = TRUE)
thr <- 0.20 * sdy   # ROPE

hdi_in_rope_pct <- function(x, thr, ci = 0.89){
  h <- bayestestR::hdi(x, ci = ci)
  lo <- h$CI_low[1]; hi <- h$CI_high[1]
  overlap <- max(0, min(hi,  thr) - max(lo, -thr))
  len     <- max(hi - lo, 1e-12)
  100 * overlap / len
}

rope_by_week <- map2_dfr(
  weeks_models$model,
  weeks_models$week_f,
  ~ {
    dr <- posterior::as_draws_df(.x)
    param_names <- grep("^b_", names(dr), value = TRUE)
    
    purrr::map_dfr(param_names, function(p) {
      vals <- as.numeric(dr[[p]])
      hd   <- bayestestR::hdi(vals, ci = 0.89)
      
      tibble(
        week_f      = .y,
        term        = p,
        mean        = mean(vals),
        hdi_low89   = hd$CI_low[1],
        hdi_high89  = hd$CI_high[1],
        hdi_in_rope = hdi_in_rope_pct(vals, thr, ci = 0.89)
      )
    })
  }
) %>%
  arrange(as.numeric(as.character(week_f)), term)

#rope_by_week
#print(rope_by_week, n = 42) # Results as shown in Table S7, Supporting Info

summary(weeks_models$model[[1]])
summary(weeks_models$model[[2]])
summary(weeks_models$model[[3]])
summary(weeks_models$model[[4]])
summary(weeks_models$model[[5]])
summary(weeks_models$model[[6]])

weeks_models2 <- weeks_models %>%
  mutate(
    preds = purrr::map2(model, data, ~ {
      dat <- .y
      log_size_seq <- seq(
        min(dat$log_size, na.rm = TRUE),
        max(dat$log_size, na.rm = TRUE),
        length.out = 100
      )
      
      temp_ref <- median(dat$temp, na.rm = TRUE)
      
      newdat <- tidyr::expand_grid(
        log_size  = log_size_seq,
        dispersal = levels(dat$dispersal)
      ) %>%
        mutate(temp = temp_ref)
      
      fitted(.x, newdata = newdat, summary = TRUE) %>%
        as_tibble() %>%
        bind_cols(newdat)
    })
  )

plot_data <- weeks_models2 %>%
  select(week_f, data) %>%
  unnest(data)

plot_preds <- weeks_models %>%
  select(week_f, preds) %>%
  unnest(preds)

plot_data <- plot_data %>%
  mutate(
    week_f = factor(
      week_f,
      levels = c("2","4","6","8","10","12"),
      labels = c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12"),
      ordered = TRUE
    )
  )

plot_preds <- plot_preds %>%
  mutate(
    week_f = factor(
      week_f,
      levels = c("2","4","6","8","10","12"),
      labels = c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12"),
      ordered = TRUE
    )
  )

week_labels <- c(
  "2"  = "Week 2",
  "4"  = "Week 4",
  "6"  = "Week 6",
  "8"  = "Week 8",
  "10" = "Week 10",  
  "12" = "Week 12"
)

disp_colors <- c("none" = "#003f5c",
                 "low"  = "#ff6361",
                 "high" = "#ffa600")

wattage_levels <- levels(phyto_size_spectra$wattage)
wattage_cols   <- brewer.pal(4, "RdBu") |> rev()  
names(wattage_cols) <- wattage_levels

# Plot size spectra by week and dispersal treatment (Figure S12a, Supporting Info)
size_spec_disp_fig <- ggplot() +
  geom_point(
    data = plot_data,
    aes(x = log_size, y = log_count, color = dispersal, shape = dispersal),
    alpha = 0.5,
    size  = 1.8,
    position = position_jitter(width = 0.02, height = 0.02)
  ) +
  geom_line(
    data = plot_preds,
    aes(x = log_size, y = Estimate, color = dispersal, linetype = dispersal),
    linewidth = 1.1
  ) +
  facet_wrap(
    ~ week_f,
    labeller = labeller(week_f = week_labels)
  ) +
  scale_linetype_manual(
    name   = "Dispersal",
    values = c("none" = "dotdash","low" = "dashed","high" = "solid"),
    labels = c("none" = "None","low" = "Low","high" = "High")
  ) +
  scale_shape_manual(
    name   = "Dispersal",
    values = c("none" = 21,"low" = 25,"high" = 24),
    labels = c("none" = "None","low" = "Low","high" = "High")
  ) +
  scale_color_manual(
    values = disp_colors,
    breaks = c("none","low","high"),
    labels = c("None","Low","High")
  ) +
  labs(
    x = expression("log cell size (ESD, " * mu * "m)"),
    y = "log abundance",
    color = "Dispersal"
  ) +
  theme_classic() +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = c("dotdash","dashed","solid"),
        shape    = c(1, 6, 2),
        fill     = c("#003f5c","#ff6361","#ffa600")
      ),
      order = 1
    ),
    linetype = "none",
    shape    = "none",
    fill     = "none"
  ) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 13, family = "Helvetica"),
    axis.title.x     = element_text(size = 13, family = "Helvetica", margin = margin(t = 11)),
    axis.title.y     = element_text(size = 13, family = "Helvetica", margin = margin(r = 11)),
    axis.text.x      = element_text(size = 10, family = "Helvetica"),
    axis.text.y      = element_text(size = 10, family = "Helvetica"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11, family = "Helvetica"),
    legend.title     = element_text(size = 12, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(0.8, "cm")
  )
size_spec_disp_fig
ggsave(size_spec_disp_fig, file="size_spectra_by_disp.pdf",width=6.1,height=4.8) #Fig. S12a




## HYPOTHESIS 2C: Thermal trait (taxonomic) composition --------------------

## PERMDISP tests on all three trophic levels, resulting in Figure 4b ##

# read in phyto count data
Phyto_data <-read.csv("phytoplankton_diversity_data.csv")
Phyto_data$tank <- as.factor(Phyto_data$tank)
metadata <- exp_data %>%
  select(sample_date, week_number,tank, wattage, metacommunity, dispersal, temp)

#Phyto_data$sample_date <- as.Date(Phyto_data$sample_date, format = "%Y-%m-%d")
# might not be necessary but make sure Phyto_data matches exp_data date format^

Phyto_combined <- Phyto_data %>%
  left_join(metadata, by = c("sample_date", "tank"))

sp_cols <- which(names(Phyto_combined) == "Ceratium") :
  which(names(Phyto_combined) == "Cryptomonas.sp.") # specify species count columns for matrix building

Phyto_combined <- Phyto_combined %>%
  tidyr::drop_na( all_of(sp_cols) ) 

Phyto_combined$week_number <- factor(Phyto_combined$week_number)

Phyto_combined$dispersal <- factor(Phyto_combined$dispersal, ordered = TRUE)
Phyto_combined$tank <- factor(Phyto_combined$tank, ordered = TRUE)

phyto_mat <- data.matrix(Phyto_combined %>% dplyr::select(Ceratium:Cryptomonas.sp.))

env_phyto <- Phyto_combined %>% # matching environment matrix for permdisp
  select(sample_date, tank, wattage, temp, week_number,dispersal,metacommunity)

env_phyto$wattage<-factor(env_phyto$wattage, levels=c("0","100","200","300"))

env_phyto <- env_phyto %>% 
  mutate(
    week = case_when(
      sample_date == as.Date("2018-07-11") ~ "Week 4",
      sample_date == as.Date("2018-08-08") ~ "Week 8",
      sample_date == as.Date("2018-09-05") ~ "Week 12",
      TRUE                                ~ NA_character_   
    ),
    week = factor(week, levels = c("Week 4", "Week 8", "Week 12"),ordered=TRUE)  
  )

##### ZOOPLANKTON #####
# Read zooplankton count data
Zoop_data <- read.csv("Zoop_data.csv")

# Get data in shape
week4zoop <- Zoop_data%>%
  filter(week==4) %>%
  dplyr::select(-c(tank,wattage,metacommunity,dispersal,week,temp_mean,date,week))

week4zoop_env <- Zoop_data %>%
  filter(week==4) %>%
  dplyr::select(c(tank,wattage,metacommunity,dispersal,week,temp_mean,date))

week8zoop <- Zoop_data %>%
  filter(week==8) %>%
  dplyr::select(-c(tank,wattage,metacommunity,dispersal,week,temp_mean,date,week))

week8zoop_env <- Zoop_data %>%
  filter(week==8) %>%
  dplyr::select(c(tank,wattage,metacommunity,dispersal,week,temp_mean,date))

week12zoop <- Zoop_data %>%
  filter(week==12) %>%
  dplyr::select(-c(tank,wattage,metacommunity,dispersal,week,temp_mean,date,week))

week12zoop_env <- Zoop_data %>%
  filter(week==12) %>%
  dplyr::select(c(tank,wattage,metacommunity,dispersal,week,temp_mean,date))

### Zoop composition PERMDISP
zoop_dist_mat4 <- vegdist(week4zoop, method = "bray")

zoop_bd4 <- betadisper(zoop_dist_mat4, type="centroid",group = week4zoop_env$metacom) 

zoop_bd_df4 <- data.frame(
  tank      = week4zoop_env$tank,
  metacom         = zoop_bd4$group, 
  dispersal= week4zoop_env$dispersal,
  DistToCentroid = zoop_bd4$distances
)

zoop_week4_dispersion <- aov(DistToCentroid ~ dispersal, data = zoop_bd_df4)
anova(zoop_week4_dispersion) 
#Response: DistToCentroid
#Df  Sum Sq   Mean Sq F value Pr(>F)
#dispersal  2 0.00300 0.0015003  0.1197 0.8875
#Residuals 45 0.56408 0.0125350   


###### now week 8 ######
zoop_dist_mat8 <- vegdist(week8zoop, method = "bray")

zoop_bd8 <- betadisper(zoop_dist_mat8, group = week8zoop_env$metacom,type="centroid")

zoop_bd_df8 <- data.frame(
  tank      = week8zoop_env$tank,
  metacom         = zoop_bd8$group, 
  dispersal= week8zoop_env$dispersal,
  DistToCentroid = zoop_bd8$distances
)

zoop_week8_dispersion <- aov(DistToCentroid ~ dispersal, data = zoop_bd_df8)
anova(zoop_week8_dispersion) 
#Response: DistToCentroid
#Df  Sum Sq   Mean Sq F value Pr(>F)
#dispersal  2 0.00728 0.0036385  0.2575 0.7741
#Residuals 45 0.63594 0.0141320  


###### now week 12 ######

zoop_dist_mat12 <- vegdist(week12zoop, method = "bray")

zoop_bd12 <- betadisper(zoop_dist_mat12, group = week12zoop_env$metacom,type="centroid")


zoop_bd_df12 <- data.frame(
  tank      = week12zoop_env$tank,
  metacom         = zoop_bd12$group, 
  dispersal= week12zoop_env$dispersal,
  DistToCentroid = zoop_bd12$distances
)

zoop_week12_dispersion <- aov(DistToCentroid ~ dispersal, data = zoop_bd_df12)
anova(zoop_week12_dispersion) 
#Response: DistToCentroid
#Df  Sum Sq  Mean Sq F value  Pr(>F)  
#dispersal  2 0.08838 0.044190  2.9437 0.06289 .
#Residuals 45 0.67554 0.015012 

## combine both
zoop_bd_df4$week  <- "Week 4"
zoop_bd_df8$week  <- "Week 8"
zoop_bd_df12$week <- "Week 12"

zoop_bd_all <- bind_rows(
  mutate(zoop_bd_df4 , week = 4),
  mutate(zoop_bd_df8 , week = 8),
  mutate(zoop_bd_df12, week = 12)
) %>% 
  relocate(week, .after = tank)  
zoop_letters_df <- zoop_bd_all %>% 
  group_by(week) %>%                    
  group_map(~ {                         
    mod <- aov(DistToCentroid ~ dispersal, data = .x)
    tuk <- TukeyHSD(mod)
    letter_vec <- multcompLetters4(mod, tuk)$dispersal$Letters
    
    tibble(
      week      = .y$week,            
      dispersal = names(letter_vec),
      letter    = as.vector(letter_vec)
    )
  }) %>% 
  bind_rows()

zoop_letters_df
zoop_bd_all$dispersal <- factor(zoop_bd_all$dispersal,levels= c("none","low","high"),ordered=TRUE)

# PHYTO COMPOSITION PERMDISP

# separate by week
week4phyto <- Phyto_combined %>%
  filter(week_number==4) %>%
  dplyr::select(Ceratium:Cryptomonas.sp.)

week4phyto_env <- Phyto_combined %>%
  filter(week_number==4) %>%
  dplyr::select(c(tank,wattage,metacommunity,dispersal,week_number,temp,sample_date))

week8phyto <- Phyto_combined %>%
  filter(week_number==8) %>%
  dplyr::select(Ceratium:Cryptomonas.sp.)

week8phyto_env <- Phyto_combined %>%
  filter(week_number==8) %>%
  dplyr::select(c(tank,wattage,metacommunity,dispersal,week_number,temp,sample_date))

week12phyto <- Phyto_combined %>%
  filter(week_number==12) %>%
  dplyr::select(Ceratium:Cryptomonas.sp.)

week12phyto_env <- Phyto_combined %>%
  filter(week_number==12) %>%
  dplyr::select(c(tank,wattage,metacommunity,dispersal,week_number,temp,sample_date))

# Run PERMANOVAs for phytoplankton taxonomic composition
phyto_perman4 <-adonis2(week4phyto ~dispersal*temp,  method = "bray", strata = week4phyto_env$metacommunity,  by="terms", data = week4phyto_env)

phyto_perman8 <-adonis2(decostand(week8phyto,method="hellinger") ~dispersal*temp,  method = "bray", strata = week8phyto_env$metacommunity,  by="terms", data = week8phyto_env)

phyto_perman12 <-adonis2(decostand(week12phyto,method="hellinger") ~dispersal*temp,  method = "bray", strata = week12phyto_env$metacommunity,  by="terms", data = week12phyto_env)


# fit PERMDISPs
phyto_dist_mat4 <- vegdist(week4phyto, method = "bray")

phyto_bd4 <- betadisper(phyto_dist_mat4, type="centroid", group = week4phyto_env$metacom) 

phyto_bd_df4 <- data.frame(
  tank      = week4phyto_env$tank,
  metacom         = phyto_bd4$group, 
  dispersal= week4phyto_env$dispersal,
  DistToCentroid = phyto_bd4$distances)

phyto_week4_dispersion <- aov(DistToCentroid ~ dispersal, data = phyto_bd_df4)
anova(phyto_week4_dispersion) 
#Response: DistToCentroid
#          Df  Sum Sq  Mean Sq F value Pr(>F)
#dispersal  2 0.02771 0.013855  0.8281 0.4439
#Residuals 42 0.70275 0.016732     


###### now week 8 ######

phyto_dist_mat8 <- vegdist(week8phyto, method = "bray")

phyto_bd8 <- betadisper(phyto_dist_mat8, group = week8phyto_env$metacom,type="centroid")

phyto_bd_df8 <- data.frame(
  tank      = week8phyto_env$tank,
  metacom         = phyto_bd8$group, 
  dispersal= week8phyto_env$dispersal,
  DistToCentroid = phyto_bd8$distances
)

phyto_week8_dispersion <- aov(DistToCentroid ~ dispersal, data = phyto_bd_df8)
anova(phyto_week8_dispersion)  
#Response: DistToCentroid
#          Df  Sum Sq   Mean Sq F value Pr(>F)
#dispersal  2 0.01135 0.0056757  0.4886 0.6168
#Residuals 44 0.51110 0.0116158    


###### now week 12 ######

phyto_dist_mat12 <- vegdist(week12phyto, method = "bray")

phyto_bd12 <- betadisper(phyto_dist_mat12, group = week12phyto_env$metacom,type="centroid")

phyto_bd_df12 <- data.frame(
  tank      = week12phyto_env$tank,
  metacom         = phyto_bd12$group, 
  dispersal= week12phyto_env$dispersal,
  DistToCentroid = phyto_bd12$distances
)

phyto_week12_dispersion <- aov(DistToCentroid ~ dispersal, data = phyto_bd_df12)
anova(phyto_week12_dispersion) 
tuk <- TukeyHSD(phyto_week12_dispersion, conf.level = 0.95)
print(tuk)

## combine 
phyto_bd_df4$week  <- "Week 4"
phyto_bd_df8$week  <- "Week 8"
phyto_bd_df12$week <- "Week 12"

phyto_bd_all <- bind_rows(
  mutate(phyto_bd_df4 , week = 4),
  mutate(phyto_bd_df8 , week = 8),
  mutate(phyto_bd_df12, week = 12)
) %>% 
  relocate(week, .after = tank)  

# get post-hoc significance letters (just for plotting purposes)
phyto_letters_df <- phyto_bd_all %>% 
  group_by(week) %>%                    
  group_map(~ {                       
    mod <- aov(DistToCentroid ~ dispersal, data = .x)
    tuk <- TukeyHSD(mod)
    letter_vec <- multcompLetters4(mod, tuk)$dispersal$Letters
    
    tibble(
      week      = .y$week,               
      dispersal = names(letter_vec),
      letter    = as.vector(letter_vec)
    )
  }) %>% 
  bind_rows()

phyto_letters_df

phyto_letters_df <- phyto_letters_df %>% 
  group_by(week) %>% 
  mutate(letter = if (all(letter == letter[1])) "a" else letter) %>% 
  ungroup()

text_pos <- phyto_bd_all %>% 
  group_by(week, dispersal) %>% 
  summarise(y = max(DistToCentroid, na.rm = TRUE) + 0.02, .groups = "drop")

phyto_letters_df <- phyto_letters_df %>% left_join(text_pos, by = c("week", "dispersal"))
phyto_letters_df <- phyto_letters_df %>%            
  left_join(
    phyto_bd_all %>%                    
      group_by(week, dispersal) %>%
      summarise(y = max(DistToCentroid, na.rm = TRUE) + 0.02,
                .groups = "drop"),
    by = c("week", "dispersal")
  )


### BACTERIAL COMPOSITION PERMDISPS ##
# read /clean bacteria data
Bacteria_dat<-read.csv("Bacterial_16S_data.csv")
Bacteria_dat$dispersal <- factor(Bacteria_dat$dispersal, ordered = TRUE)
Bacteria_dat$tank <- as.factor(Bacteria_dat$tank)

# Get data in shape
week4bacteria <- Bacteria_dat%>%
  filter(week==4) %>%
  dplyr::select(ASV1:ASV1413)

week4bacteria_env <- Bacteria_dat %>%
  filter(week==4) %>%
  dplyr::select(c(tank,wattage,metacommunity,dispersal,week,temp_mean))

week12bacteria <- Bacteria_dat %>%
  filter(week==12) %>%
  dplyr::select(ASV1:ASV1413)

week12bacteria_env <- Bacteria_dat %>%
  filter(week==12) %>%
  dplyr::select(c(tank,wattage,metacommunity,dispersal,week,temp_mean))

bacteria_dist_mat4 <- vegdist(week4bacteria, method = "bray")

bacteria_bd4 <- betadisper(bacteria_dist_mat4, type="centroid",group = week4bacteria_env$metacom) 

bacteria_bd_df4 <- data.frame(
  tank      = week4bacteria_env$tank,
  metacom         = bacteria_bd4$group, 
  dispersal= week4bacteria_env$dispersal,
  DistToCentroid = bacteria_bd4$distances
)

bacteria_week4_dispersion <- aov(DistToCentroid ~ dispersal, data = bacteria_bd_df4)
anova(bacteria_week4_dispersion) # significant week 4

#Response: DistToCentroid
#Df   Sum Sq  Mean Sq F value  Pr(>F)  
#dispersal  2 0.025403 0.012701  3.4543 0.04063 *
#  Residuals 43 0.158112 0.003677         


###### now week 12 ######

bacteria_dist_mat12 <- vegdist(week12bacteria, method = "bray")

bacteria_bd12 <- betadisper(bacteria_dist_mat12, group = week12bacteria_env$metacom,type="centroid")

# Make a data frame from the betadisper object
bacteria_bd_df12 <- data.frame(
  tank      = week12bacteria_env$tank,
  metacom         = bacteria_bd12$group, 
  dispersal= week12bacteria_env$dispersal,
  DistToCentroid = bacteria_bd12$distances
)

bacteria_week12_dispersion <- aov(DistToCentroid ~ dispersal, data = bacteria_bd_df12)
anova(bacteria_week12_dispersion) 
tuk <- TukeyHSD(bacteria_week12_dispersion, conf.level = 0.95)
print(tuk)
#Response: DistToCentroid
#Df   Sum Sq  Mean Sq F value   Pr(>F)   
#dispersal  2 0.049674 0.024837  8.2488 0.001029 **
#  Residuals 39 0.117428 0.003011      

## combine 

bacteria_bd_all <- bind_rows(
  mutate(bacteria_bd_df4 , week = 4),
  mutate(bacteria_bd_df12, week = 12)
) %>% 
  relocate(week, .after = tank)  
bacteria_bd_all$dispersal<-factor(bacteria_bd_all$dispersal,levels=c("none","low","high"))

bacteria_letters_df <- bacteria_bd_all %>% 
  group_by(week) %>%                     
  group_map(~ {                         
    mod <- aov(DistToCentroid ~ dispersal, data = .x)
    tuk <- TukeyHSD(mod)
    letter_vec <- multcompLetters4(mod, tuk)$dispersal$Letters
    
    tibble(
      week      = .y$week,              
      dispersal = names(letter_vec),
      letter    = as.vector(letter_vec)
    )
  }) %>% 
  bind_rows()

bacteria_letters_df <- bacteria_letters_df %>% 
  group_by(week) %>% 
  mutate(letter = if (all(letter == letter[1])) "a" else letter) %>% 
  ungroup()

text_pos <- bacteria_bd_all %>% 
  group_by(week, dispersal) %>% 
  summarise(y = max(DistToCentroid, na.rm = TRUE) + 0.02, .groups = "drop")

bacteria_letters_df <- bacteria_letters_df %>% left_join(text_pos, by = c("week", "dispersal"))
bacteria_letters_df <- bacteria_letters_df %>%           
  left_join(
    bacteria_bd_all %>%                   
      group_by(week, dispersal) %>%
      summarise(y = max(DistToCentroid, na.rm = TRUE) + 0.02,
                .groups = "drop"),
    by = c("week", "dispersal")
  )

### Now combine all for manuscript Figure 4b

disp_pal <- c("#003f5c", "#ff6361", "#ffa600")

LEG_TXT  <- 13.5
LEG_TTL  <- 14.5
XTICK    <- 13.5
XTITLE   <- 14.5
TITLE_SZ <- 14
YLAB     <- "Distance to metacommunity centroid"

common_fill <- scale_fill_manual(values = disp_pal, name = "Dispersal", limits = c("none","low","high"))

disp_colors <- c("none" = "#003f5c",
                 "low"  = "#ff6361",
                 "high" = "#ffa600")
title_theme <- theme(
  plot.title.position = "plot",
  plot.title = element_text(size = TITLE_SZ, hjust = 0.5, margin = margin(l = 8))
)

summarise_dispersion <- function(df) {
  df %>%
    mutate(
      week_f = factor(week, levels = c(4, 8, 12),
                      labels = c("4", "8", "12"))
    ) %>%
    group_by(week_f, dispersal) %>%
    summarise(
      mean_dist = mean(DistToCentroid, na.rm = TRUE),
      se_dist   = sd(DistToCentroid, na.rm = TRUE) /
        sqrt(sum(!is.na(DistToCentroid))),
      .groups   = "drop"
    )
}

phyto_sum    <- summarise_dispersion(phyto_bd_all)
bacteria_sum <- summarise_dispersion(bacteria_bd_all)
zoop_sum     <- summarise_dispersion(zoop_bd_all)

p_phyto_pts <- ggplot(
  phyto_sum,
  aes(x = week_f, y = mean_dist, color = dispersal, group = dispersal)
) +geom_point(position = position_dodge(width = 0.5), size = 3,aes(shape=dispersal)) +
  geom_errorbar(
    aes(ymin = mean_dist - se_dist,
        ymax = mean_dist + se_dist),
    width     = 0.3,
    position  = position_dodge(width = 0.5),
    linewidth = 0.6
  ) +
  scale_color_manual(values=disp_colors) + # if this sets scale_color / legends etc.
  scale_y_continuous(
    limits = c(0.23, 0.5),
    breaks = c(0.3, 0.4),
    expand = c(0, 0)
  ) +
  scale_x_discrete(name = "", labels = c("4","8","12")) +
  labs(y = "") +
  theme_classic(base_size = 13) +
  title_theme +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = LEG_TXT),
    legend.title    = element_text(size = LEG_TTL),
    panel.border    = element_rect(colour = "black", fill = NA),
    axis.text.x     = element_text(size = XTICK),
    axis.text.y     = element_text(size = 12),
    axis.title.y    = element_text(size = 14, margin = margin(r = 10)),
    plot.margin     = margin(22, 10, 6, 6)
  )

p_bacteria_pts <- ggplot(
  bacteria_sum,
  aes(x = week_f, y = mean_dist, color = dispersal, group = dispersal)
) +
  geom_point(position = position_dodge(width = 0.5), size = 3,aes(shape=dispersal)) +
  geom_errorbar(
    aes(ymin = mean_dist - se_dist,
        ymax = mean_dist + se_dist),
    width     = 0.3,
    position  = position_dodge(width = 0.5),
    linewidth = 0.6
  ) +
  scale_color_manual(values=disp_colors) +
  scale_y_continuous(
    limits = c(0.28 , 0.43),
    breaks = c(0.3, 0.4),
    expand = c(0, 0)
  ) +
  scale_x_discrete(name = "", labels = c("4","8","12")) +
  labs(y = "") +
  theme_classic(base_size = 13) +
  title_theme +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = LEG_TXT),
    legend.title    = element_text(size = LEG_TTL),
    panel.border    = element_rect(colour = "black", fill = NA),
    axis.text.x     = element_text(size = XTICK),
    axis.text.y     = element_text(size = 12),
    axis.title.x    = element_text(size = XTITLE),
    plot.margin     = margin(6, 10, 5, 6)
  )


p_zoop_pts <- ggplot(
  zoop_sum,
  aes(x = week_f, y = mean_dist, color = dispersal, group = dispersal)
) +
  geom_point(position = position_dodge(width = 0.5), size = 3,aes(shape=dispersal)) +
  geom_errorbar(
    aes(ymin = mean_dist - se_dist,
        ymax = mean_dist + se_dist),
    width     = 0.3,
    position  = position_dodge(width = 0.5),
    linewidth = 0.6
  ) +
  scale_color_manual(values=disp_colors) +
  scale_y_continuous(
    limits = c(0.1, 0.5),
    expand = c(0, 0)
  ) +
  scale_x_discrete(name = "", labels = c("4","8","12"))+
  labs(y = "") +
  theme_classic(base_size = 13) +
  title_theme +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = LEG_TXT),
    legend.title    = element_text(size = LEG_TTL),
    panel.border    = element_rect(colour = "black", fill = NA),
    axis.text.x     = element_text(size = XTICK),
    axis.text.y     = element_text(size = 12),
    axis.title.x    = element_text(size = XTITLE, vjust = 1),
    plot.margin     = margin(6, 10, 10, 10)
  )

# Putting Fig 4b together:
final_permdisp_plot <-
  (p_bacteria_pts | p_phyto_pts | p_zoop_pts) +
  plot_layout(widths = c(1, 1,1), guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.margin = margin(15,15,15,15)  
  ) &
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


#ggsave(final_permdisp_plot,"permdisp.pdf",width=8.4,height=7.8)



####
### Hypothesis 2c continued: 
### Estimating and visualizing thermal phenotype composition (Fig. 5 and S13) ###

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(rTPC)

exp_data <- read.csv("metacom_exp_data_final.csv")
therm_traits<-read.csv("trait_database.csv")
Phyto_data <- read.csv("phytoplankton_diversity_data.csv")

disp_colors <- c(
  "none" = "#003f5c",
  "low"  = "#ff6361",
  "high" = "#ffa600"
) # for plotting later

## Data prep: taking mean trait value for each taxon if there are more than one values

traits_summ<-therm_traits %>%
  group_by(genus) %>%
  summarize(across(c(w, a, b, z, T_opt, max_growth_at_Topt,tmin,tmax), ~ mean(.x, na.rm = TRUE)))

traits_summ_scaled <- traits_summ %>%
  mutate(
    genus_clean = genus %>%
      str_replace_all("\\.", " ") %>%
      str_squish(),
    breadth = tmax -
      tmin
  )


# Plotting TPCs for phyto taxa found in experiment (Fig. S13a) 

therm_traits_tpc <- traits_summ %>%
  mutate(
    breadth = ifelse(!is.na(tmin) & !is.na(tmax), tmax - tmin, NA_real_)
  )
fallback_breadth <- median(therm_traits_tpc$breadth, na.rm = TRUE)
fallback_height  <- median(therm_traits_tpc$max_growth_at_Topt, na.rm = TRUE)

# temperature range for plotting
T_min <- min(therm_traits_tpc$tmin, therm_traits_tpc$T_opt, na.rm = TRUE) - 2
T_max <- max(therm_traits_tpc$tmax, therm_traits_tpc$T_opt, na.rm = TRUE) + 2
T_seq <- seq(T_min, T_max, length.out = 200)

# Thomas 2012 TPC model helper from Daniel Padfield's rTPC package
make_tpc <- function(row) {
  row <- as.list(row)
  
  if (!is.na(row$T_opt) && !is.na(row$max_growth_at_Topt) && !is.na(row$breadth)) {
    a    <- row$max_growth_at_Topt
    b    <- 0
    cpar <- row$breadth
    topt <- row$T_opt
    
    rate <- rTPC::thomas_2012(temp = T_seq, a = a, b = b, c = cpar, topt = topt)
    
    tibble(
      T       = T_seq,
      rate    = rate,
      species = row$species,
      genus   = row$genus,
      model   = "thomas_full"
    )
    
  } else if (!is.na(row$T_opt) && !is.na(row$max_growth_at_Topt)) {
    
    a    <- row$max_growth_at_Topt
    b    <- 0
    cpar <- ifelse(is.na(row$breadth), fallback_breadth, row$breadth)
    topt <- row$T_opt
    
    rate <- rTPC::thomas_2012(temp = T_seq, a = a, b = b, c = cpar, topt = topt)
    
    tibble(
      T       = T_seq,
      rate    = rate,
      species = row$species,
      genus   = row$genus,
      model   = "thomas_partial"
    )
    
  } else if (!is.na(row$T_opt)) {
    
    sigma <- 2  
    rate  <- dnorm(T_seq, mean = row$T_opt, sd = sigma) *
      ifelse(is.na(row$max_growth_at_Topt), fallback_height, row$max_growth_at_Topt)
    
    tibble(
      T       = T_seq,
      rate    = rate,
      species = row$species,
      genus   = row$genus,
      model   = "gaussian_fallback"
    )
    
    
  } else {
    tibble()
  }
}

# generate TPCs for each species (Fig 13a)
tpcs_thomas <- therm_traits_tpc %>%
  group_by(genus) %>%
  slice(1) %>%
  ungroup() %>%
  split(.$genus) %>%
  map_dfr(~ make_tpc(.x[1, ]))

# plot
species_tpcs <- ggplot(tpcs_thomas, aes(x = T, y = rate, group = genus, color = genus)) +
  geom_line(size=0.8,alpha = 0.7) +
  theme_bw() +
  labs(
    x = "Temperature (°C)",
    y = "Population growth rate (modelled)"
  )+xlim(0,40)+ylim(0,1.75)+theme_classic()+
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
    axis.text.x      = element_text(size = 11, family = "Helvetica"),
    axis.text.y      = element_text(size = 11, family = "Helvetica"),
    legend.text      = element_text(size = 11, family = "Helvetica"),
    legend.title     = element_text(size = 12, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(0.8, "cm")
  )

#ggsave(species_tpcs,file="species_tpcs.pdf",width=5.8,height=3.4) # as in Fig S13a

### Looking for correlation between Topt and max growth rates (hotter is better vs compensation) --> Fig. S13b in Supp Info

df <- therm_traits %>%
  select(T_opt, max_growth_at_Topt,w) %>%
  filter(!is.na(T_opt), !is.na(max_growth_at_Topt),!is.na(w))

lm_lin <- lm(max_growth_at_Topt ~ T_opt, data = df)
summary(lm_lin)

lm_exp <- lm(log(max_growth_at_Topt) ~ T_opt, data = df)
summary(lm_exp)

#Compare model fit
AIC(lm_lin, lm_exp)

# plot 
max_v_topt<- ggplot(df, aes(x = T_opt, y = log(max_growth_at_Topt))) +
  geom_point(size = 2) + geom_smooth(method="lm",se=TRUE)+
  labs(
    y = "log max growth rate at Topt",
    x = "Topt (°C)"
  )+theme_classic()+
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 10),
    axis.title       = element_text(size = 12),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    legend.title     = element_text(size = 12)
  )

#ggsave(max_v_topt, file="max_v_topt.pdf",width=4.2,height=3.4) # as in Fig S13b

# Calculating CTPI (Community Thermal Phenotype Index)

# standardize trait ranges
Topt_range   <- range(traits_summ_scaled$T_opt,   na.rm = TRUE)
breadth_range <- range(traits_summ_scaled$breadth, na.rm = TRUE)

traits_summ_scaled <- traits_summ_scaled %>%
  mutate(
    Topt_s = if (diff(Topt_range) == 0) 0.5 else
      (T_opt   - Topt_range[1])   / diff(Topt_range),
    breadth_s = if (diff(breadth_range) == 0) 0.5 else
      (breadth - breadth_range[1]) / diff(breadth_range),
    trait_score = Topt_s * breadth_s
  ) %>%
  select(genus_clean, T_opt, tmin, tmax, breadth,
         Topt_s, breadth_s, trait_score)

Phyto_data <- Phyto_data %>%
  dplyr::rename(
    Oocystis.large = Large.oocystis,
    Oocystis.small = Small.oocystis
  )

phyto_long <- Phyto_data %>%
  pivot_longer(
    cols = -c(sample_date, tank),
    names_to = "species_raw",
    values_to = "count"
  ) %>%
  filter(!is.na(count)) %>%
  mutate(
    genus_clean = species_raw %>%
      sub("\\..*$", "", .) %>%       # "Acutodesmus.obliquus" -> "Acutodesmus"
      str_replace_all("\\.", " ") %>%
      str_squish()
  )

# Join back with metadata
comm_traits <- phyto_long %>%
  left_join(traits_summ_scaled, by = "genus_clean") %>%
  filter(!is.na(trait_score))  

thermal_index_by_tank_date <- comm_traits %>%
  group_by(tank, sample_date) %>%
  mutate(
    rel_abund = count / sum(count, na.rm = TRUE)
  ) %>%
  summarise(
    CTPI_raw = sum(rel_abund * trait_score, na.rm = TRUE),
    n_genera_with_traits = n_distinct(genus_clean),
    .groups = "drop"
  )

### Rescale CTPI_raw to 0–1 across all tank × dates 
CTPI_range <- range(thermal_index_by_tank_date$CTPI_raw, na.rm = TRUE)

thermal_index_by_tank_date <- thermal_index_by_tank_date %>%
  mutate(
    CTPI_scaled = if (diff(CTPI_range) == 0) 0.5 else
      (CTPI_raw - CTPI_range[1]) / diff(CTPI_range)
  )

thermal_index_by_tank_date$tank <- as.factor(thermal_index_by_tank_date$tank)

exp_data <- exp_data %>%
  left_join(
    thermal_index_by_tank_date %>% select(tank, sample_date, CTPI_scaled),
    by = c("tank", "sample_date")
  )


ctpi_sum <- exp_data %>%
  filter(week_number %in% c(4, 8, 12)) %>%
  group_by(week_number, dispersal) %>%
  summarise(
    mean_CTPI = mean(CTPI_scaled, na.rm = TRUE),
    se_CTPI   = sd(CTPI_scaled, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

# Plotting CTPI by dispersal treatment over weeks 4, 8, 12
pd  <- position_dodge(width = 0.3)
pjd <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.3)

CTPI_summary <- ggplot() +
  geom_point(
    data = exp_data %>% filter(week_number %in% c(4, 8, 12)),
    aes(
      x     = factor(week_number),
      y     = CTPI_scaled,
      color = dispersal
    ),
    position = pjd,
    alpha    = 0.2,
    size     = 1.7
  ) +
  geom_errorbar(
    data = ctpi_sum,
    aes(
      x    = factor(week_number),
      y    = mean_CTPI,
      ymin = mean_CTPI - se_CTPI,
      ymax = mean_CTPI + se_CTPI,
      color = dispersal
    ),
    position  = pd,
    width     = 0.25,
    linewidth = 0.7
  ) +
  geom_point(
    data = ctpi_sum,
    aes(
      x     = factor(week_number),
      y     = mean_CTPI,
      color = dispersal,
      shape = dispersal
    ),
    position = pd,
    size     = 3
  ) +
  
  scale_color_manual(values = disp_colors) +
  labs(
    y = "Phytoplankton Community Thermal Phenotype Index",
    x = "Week"
  ) +
  theme_classic() +
  theme(
    legend.position  = "none",
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
    axis.text.x      = element_text(size = 11, family = "Helvetica"),
    axis.text.y      = element_text(size = 11, family = "Helvetica"),
    legend.text      = element_text(size = 11, family = "Helvetica"),
    legend.title     = element_text(size = 12, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(0.8, "cm")
  )

# Note this figure isn't in main text or supplement
#ggsave(CTPI_summary,file="CTPI_summary.pdf",width=4.7,height=4.2)

## #Tukey test to see if CTPI varies with dispersal within each week:

exp_data$week_number <- as.factor(exp_data$week_number)

test_CTPI_per_week <- function(week){
  df <- exp_data %>% filter(week_number == week)
  
  model <- lm(CTPI_scaled ~ dispersal, data = df)
  anova_table <- anova(model)
  
  posthoc <- emmeans(model, pairwise ~ dispersal)$contrasts
  
  list(
    week = week,
    anova = anova_table,
    posthoc = posthoc
  )
}

# Run for each week
results_week4  <- test_CTPI_per_week("4")
results_week8  <- test_CTPI_per_week("8")
results_week12 <- test_CTPI_per_week("12")

results_week4
results_week8
results_week12

### Calculating Community Thermal Match Index (CTMI) 

traits_summ_clean <- traits_summ %>%
  mutate(
    genus_clean = genus %>%
      str_replace_all("\\.", " ") %>%  # just in case
      str_squish()
  ) %>%
  select(genus_clean, T_opt)

### reshape
phyto_long <- Phyto_data %>%
  pivot_longer(
    cols = -c(sample_date, tank),
    names_to = "species_raw",
    values_to = "count"
  ) %>%
  filter(!is.na(count)) %>%
  mutate(
    genus_clean = species_raw %>% # cleaning taxon names
      sub("\\..*$", "", .) %>%          
      str_replace_all("\\.", " ") %>%   
      str_squish()
  )

### Clean exp_data date format + keep temperature (unique tank-date)
exp_clean <- exp_data %>%
  select(tank, sample_date, temp) %>%
  distinct()

### rejoin w community data

phyto_long$tank <- as.factor(phyto_long$tank)

comm_traits <- phyto_long %>%
  left_join(traits_summ_clean, by = "genus_clean") %>%
  left_join(exp_data, by = c("tank", "sample_date")) %>%
  filter(!is.na(T_opt), !is.na(temp))

###  Compute CTMI (Community Thermal Mismatch Index) per tank × sample_date
CTMI_by_tank_date <- comm_traits %>%
  group_by(tank, sample_date) %>%
  mutate(
    rel_abund = count / sum(count, na.rm = TRUE)   # relative abundance within tank&date
  ) %>%
  summarise(
    CTMI_raw = sum(rel_abund * abs(T_opt - temp), na.rm = TRUE),  # abundance-weighted |Topt - T|
    n_genera_with_Topt = n_distinct(genus_clean),
    .groups = "drop"
  )

### Scale CTMI_raw to 0–1
range_vals <- range(CTMI_by_tank_date$CTMI_raw, na.rm = TRUE)

if (diff(range_vals) == 0) {
  CTMI_by_tank_date <- CTMI_by_tank_date %>%
    mutate(CTMI_scaled = 1)
} else {
  CTMI_by_tank_date <- CTMI_by_tank_date %>%
    mutate(
      CTMI_scaled = (range_vals[2] - CTMI_raw) / diff(range_vals)
    )
}

# rejoin with expdata
exp_data <- exp_data %>%
  left_join(
    CTMI_by_tank_date %>% select(tank, sample_date, CTMI_scaled),
    by = c("tank", "sample_date")
  )

##### CTMI by wattage and dispersal over time (Fig. 5b) ####

# datawrangling
ctmi_sum <- exp_data %>%
  filter(week_number %in% c(4, 8, 12)) %>%
  group_by(week_number, wattage, dispersal) %>%
  summarise(
    mean_CTMI = mean(CTMI_scaled, na.rm = TRUE),
    se_CTMI   = sd(CTMI_scaled, na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  )

exp_data <- exp_data %>%
  mutate(
    week_f = factor(
      week_number,
      levels = c(4, 8, 12),
      labels = c("Week 4", "Week 8", "Week 12")
    )
  )

ctmi_sum <- ctmi_sum %>%
  mutate(
    week_f = factor(
      week_number,
      levels = c(4, 8, 12),
      labels = c("Week 4", "Week 8", "Week 12")
    )
  )
#plot Figure 5b
CTMI_fig5b <- ggplot() +
  geom_jitter(
    data = exp_data %>% filter(week_number %in% c(4,8,12)),
    aes(x = wattage, y = CTMI_scaled, color = dispersal),
    width  = 0.15,
    height = 0,
    alpha  = 0.4,
    size   = 1.6
  ) +
  
  geom_errorbar(
    data = ctmi_sum,
    aes(
      x    = wattage,
      y    = mean_CTMI,
      ymin = mean_CTMI - se_CTMI,
      ymax = mean_CTMI + se_CTMI,
      color = dispersal,
    ),
    width = 0.3,
    linewidth = 0.7,
    position = position_dodge(width = 0.5)
  ) +
  geom_point(
    data = ctmi_sum,
    aes(x = wattage, y = mean_CTMI, color = dispersal, shape = dispersal),
    size = 3,
    position = position_dodge(width = 0.5),
  ) +
  geom_smooth(
    data = exp_data %>% filter(week_number %in% c(4,8,12)),
    aes(x = wattage, y = CTMI_scaled, color = dispersal),
    method = "gam",
    se = FALSE,
    linewidth = 1
  ) +
  
  facet_grid(~ week_f, scales = "free_y") +
  scale_color_manual(values = disp_colors) +
  theme_classic() +
  labs(
    x = "Wattage",
    y = "Community Thermal Match Index",
    color = "Dispersal",
    shape = "Dispersal"
  ) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 16, family = "Helvetica"),
    axis.title.x     = element_text(size = 16, family = "Helvetica",
                                    margin = margin(t = 11)),
    axis.title.y     = element_text(size = 16, family = "Helvetica",
                                    margin = margin(r = 11)),
    axis.text.x      = element_text(size = 11, family = "Helvetica"),
    axis.text.y      = element_text(size = 11, family = "Helvetica"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12, family = "Helvetica"),
    legend.title     = element_text(size = 12, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(0.8, "cm")
  )





