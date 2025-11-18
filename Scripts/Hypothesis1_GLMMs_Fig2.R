library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(brms)
library(emmeans)
library(tidybayes)
library(stringr)
library(forcats)
library(tibble)

exp_data<-read.csv("metacom_exp_data_final.csv")

exp_data$metacommunity <- factor(exp_data$metacommunity)
exp_data$wattage <- factor(exp_data$wattage, levels=c("0","100","200","300"))
exp_data$dispersal <- factor(exp_data$dispersal, levels=c("none","low","high"))
exp_data$tank <- factor(exp_data$tank)
exp_data$sample_date <- as.character(exp_data$sample_date)
#unique(exp_data$week_number)

dat0 <- exp_data %>%
  mutate(
    week_number   = factor(
      week_number,
      levels = c("2","4","6","8","10","12"),
      labels = c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12"),
      ordered = TRUE
    ),
    dispersal     = factor(dispersal, levels = c("none","low","high")),
    metacommunity = factor(metacommunity),
    tank          = factor(tank)
  )

WLVL <- c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12")

dat0 <- dat0 %>%
  mutate(
    ln_fluoro_ug_C_L = scale(ln_fluoro_ug_C_L)
  )

priors <- c(
  prior(normal(1.5,5),class = "Intercept"),
  prior(normal(0.65, 0.28),class = "b", coef = "inv_kTc"),
  prior(normal(0, 0.5), class = "b", coef = "inv_kTc:dispersallow"),
  prior(normal(0, 0.5), class = "b", coef = "inv_kTc:dispersalhigh"),
  prior(normal(0.75, 0.25), class = "b", coef = "ln_fluoro_ug_C_L"), # around 0.75 scaling coeff
  prior(exponential(1),   class = "sigma"),
  prior(exponential(1),   class = "sd")
)


# Fit the model for each week sampled (week_number)
# GLMM is Log(GPP) = temp + tempXdisp + ln biomass + metacommunity random effect
# where temp = 1/kTc-1/kT (standardized arrhenius temperature)

fit_one_week <- function(d) {
  dat_d <- dat0 %>%
    filter(week_number == d) %>%
    tidyr::drop_na(ln_GPP_umol_L_h, inv_kTc, ln_fluoro_ug_C_L, dispersal, metacommunity, tank)
  
  stopifnot(nrow(dat_d) >= 8)
  
  brm(
    ln_GPP_umol_L_h ~ inv_kTc + inv_kTc:dispersal + ln_fluoro_ug_C_L + (1 |metacommunity), 
    data      = dat_d,
    family    = gaussian(),
    prior     = priors,
    chains    = 4, cores = 4, iter = 5000,
    backend   = "cmdstanr",
    control   = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    seed      = 123
  )
}

# Fit all weeks (names are strings matching WLVL)
weeks <- levels(dat0$week_number)
fits  <- setNames(map(weeks, fit_one_week), weeks)

disp_levels <- c("none","low","high")
xr <- range(dat0$inv_kTc, na.rm = TRUE)
yr <- c(2.2, 3.2)

x_pad <- 0.006 * diff(xr)
y_pad <- 0.05 * diff(yr)
y_gap <- 0.2 * diff(yr)

lab_pos <- tibble(
  dispersal = factor(disp_levels, levels = disp_levels),
  x = xr[1] + x_pad-0.15,
  y = (yr[2] - y_pad) + 0.3 - (seq_along(disp_levels) - 1) * y_gap
)
.temp_slopes_from_fit <- function(fit) {
  d  <- posterior::as_draws_df(fit)
  bT   <- d$`b_inv_kTc`
  bTL  <- if ("b_inv_kTc:dispersallow"  %in% names(d)) d$`b_inv_kTc:dispersallow`  else 0
  bTH  <- if ("b_inv_kTc:dispersalhigh" %in% names(d)) d$`b_inv_kTc:dispersalhigh` else 0
  
  tibble(
    dispersal  = factor(c("none","low","high"), levels = c("none","low","high")),
    slope_mean = c(mean(bT), mean(bT + bTL), mean(bT + bTH))
  )
}

slopes_by_week <- map_dfr(names(fits), function(ww) {
  .temp_slopes_from_fit(fits[[ww]]) %>%
    mutate(week_number = factor(ww, levels = WLVL, ordered = TRUE))
}) %>%
  select(week_number, dispersal, slope_mean) %>%
  arrange(week_number, dispersal)


slopes_lab <- slopes_by_week %>%
  left_join(lab_pos, by = "dispersal") %>%
  mutate(label = paste0(dispersal, ": ", sprintf("%.2f", slope_mean), " eV"))

## fitting regression lines
pred_lines <- map_dfr(names(fits), function(ww) {
  fit  <- fits[[ww]]
  dat_w <- dat0 %>% filter(week_number == ww) %>%
    drop_na(inv_kTc, ln_fluoro_ug_C_L, dispersal)
  
  xr_w <- range(dat_w$inv_kTc, na.rm = TRUE)
  inv_seq <- seq(xr_w[1] - 0.5, xr_w[2] + 0.5, length.out = 160)  
  lnC_fix <- mean(dat_w$ln_fluoro_ug_C_L, na.rm = TRUE)
  
  newd <- expand.grid(inv_kTc = inv_seq, dispersal = levels(dat_w$dispersal)) %>%
    tibble::as_tibble() %>%
    mutate(
      ln_fluoro_ug_C_L = lnC_fix,
      week_number      = factor(ww, levels = WLVL, ordered = TRUE),
      metacommunity    = dat_w$metacommunity[1],
      tank             = dat_w$tank[1],
      dispersal_plot   = as.character(dispersal) 
    )
  
  mu <- posterior_epred(fit, newdata = newd, re_formula = NA)
  
  newd %>%
    bind_cols(as_tibble(t(mu))) %>%
    pivot_longer(starts_with("V"), values_to = "draw", names_to = "draw_id") %>%
    group_by(week_number, dispersal, dispersal_plot, inv_kTc) %>%
    summarise(
      y_mean = mean(draw),
      y_lo   = quantile(draw, 0.025),
      y_hi   = quantile(draw, 0.975),
      .groups = "drop"
    )
})
pred_lines_plot <- pred_lines

## Generate plot
dat0 <- dat0 %>%
  mutate(dispersal_plot = factor(dispersal, levels = c("none","low","high")))

fig4 <- ggplot() +
  geom_point(
    data = dat0,
    aes(inv_kTc, ln_GPP_umol_L_h,
        color  = dispersal_plot,
        shape = dispersal_plot),
    size  = 2.5,          
    alpha = 0.8     
  ) +
  geom_ribbon(
    data = pred_lines_plot,
    aes(inv_kTc, ymin = y_lo, ymax = y_hi, fill = dispersal_plot),
    alpha = 0.07
  ) +
  geom_line(
    data = pred_lines_plot,
    aes(inv_kTc, y = y_mean, color = dispersal_plot, linetype = dispersal_plot),
    linewidth = 1.3
  ) +
  geom_text(
    data = slopes_lab,
    aes(x = x, y = y, label = label),
    hjust = 0, size = 4.3, color = "black"
  ) +
  facet_wrap(~ week_number, nrow = 1, drop = TRUE) +
  guides(alpha = "none", shape = "none", fill = "none") +
  scale_color_manual(
    name   = "Dispersal",
    values = c("none"="#003f5c","low"="#ff6361","high"="#ffa600"),
    labels = c("none"="None","low"="Low","high"="High")
  ) +
  scale_linetype_manual(
    name   = "Dispersal",
    values = c("none"="dotdash","low"="dashed","high"="solid"),
    labels = c("none"="None","low"="Low","high"="High")
  ) +
  scale_shape_manual(
    name   = "Dispersal",
    values = c("none" = 21,  
               "low"  = 25,  
               "high" = 24), 
    labels = c("none" = "None",
               "low"  = "Low",
               "high" = "High")
  )+
  scale_fill_manual(  
    values = c("none"="#003f5c","low"="#ff6361","high"="#ffa600"),
    labels = c("none"="None","low"="Low","high"="High"),
    guide  = "none"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = c("dotdash","dashed","solid"),
        shape    = c(1,6,2),
        fill     = c("#003f5c", "#ff6361", "#ffa600")  
      ),
      order = 1
    ),
    linetype = "none",
    shape    = "none",
    fill     = "none"
  ) +
  labs(
    x = expression("Standardized temperature, " * 1 / italic(k)*italic(T)[c] - 1 / italic(k)*italic(T)*" (1/eV)"),
    y = expression(ln~GPP~"(µmol L"^"-1"~"h"^"-1"*")"),
    color = "Dispersal", linetype = "Dispersal", shape = "Dispersal"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "white", color = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 16, family = "Helvetica"),
    axis.title.x     = element_text(size = 16, family = "Helvetica", margin = margin(t = 11)),
    axis.title.y     = element_text(size = 16, family = "Helvetica", margin = margin(r = 11)),
    axis.text.x      = element_text(size = 15, family = "Helvetica"),
    axis.text.y      = element_text(size = 15, family = "Helvetica"),
    legend.text      = element_text(size = 13, family = "Helvetica"),
    legend.title     = element_text(size = 14, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(1, "cm")
  )
fig4

dat0wk10 <- dat0 %>% filter(week_number=="Week 10")
pred_lines_plot_wk10 <- pred_lines_plot %>% filter(week_number==10)

ggsave(fig4a, file="gpp_v_temp.pdf",width=14.5,height = 3.8)

### Supplementary version with biomass on X axis instead of temperature
dat_for_plot <- dat0 %>%
  tidyr::drop_na(ln_GPP_umol_L_h, inv_kTc, ln_fluoro_ug_C_L, dispersal, week_number) %>%
  mutate(
    dispersal_plot = factor(dispersal, levels = c("none","low","high"),
                            labels = c("none","low","high"))
  )


slopes_biomass_by_week <- purrr::map_dfr(names(fits), function(ww) {
  fit <- fits[[ww]]
  em  <- emmeans::emtrends(fit, ~ 1, var = "ln_fluoro_ug_C_L", re_formula = NA)
  edf <- as.data.frame(summary(em))
  
  est_col   <- if ("emtrends" %in% names(edf)) "emtrends" else grep("trend$", names(edf), value = TRUE)[1]
  lower_col <- if ("lower.HPD" %in% names(edf)) "lower.HPD" else "lower.CL"
  upper_col <- if ("upper.HPD" %in% names(edf)) "upper.HPD" else "upper.CL"
  
  tibble::tibble(
    week_number = factor(ww, levels = WLVL, ordered = TRUE),
    slope_mean  = edf[[est_col]],
    slope_lower = edf[[lower_col]],
    slope_upper = edf[[upper_col]]
  )
}) %>% arrange(week_number)

# flag if 95% CI excludes 0
biomass_sig_flag <- !(slopes_biomass_by_week$slope_lower <= 0 &
                        slopes_biomass_by_week$slope_upper >= 0)

slopes_lab_biomass <- slopes_biomass_by_week %>%
  mutate(
    lbl_core = paste0("biomass: ", sprintf("%.2f", slope_mean),
                      ifelse(biomass_sig_flag, " *", "")) 
  )

ranges_by_week <- dat_for_plot %>%
  dplyr::group_by(week_number) %>%
  dplyr::summarise(
    x_min = min(ln_fluoro_ug_C_L, na.rm = TRUE),
    x_max = max(ln_fluoro_ug_C_L, na.rm = TRUE),
    y_min = min(ln_GPP_umol_L_h,  na.rm = TRUE),
    y_max = max(ln_GPP_umol_L_h,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    x_lab = x_min + 0.02 * (x_max - x_min),
    y_top = y_max - 0.02 * (y_max - y_min)
  )

disp_levels <- c("none","low","high")
y_gap_frac  <- 0.10  

biomass_glmm <- ggplot() +
  geom_point(data = dat0,
             aes(x = ln_fluoro_ug_C_L, y = ln_GPP_umol_L_h,
                 color = wattage, shape = dispersal_plot),
             size = 2.5, alpha = 0.4) +
  geom_ribbon(data = pred_lines_biomass,
              aes(x = ln_fluoro_ug_C_L, ymin = y_lo, ymax = y_hi),
              alpha = 0.08) +
  geom_line(data = pred_lines_biomass,
            aes(x = ln_fluoro_ug_C_L, y = y_mean,
                alpha=0.7, linetype = dispersal_plot),
            linewidth = 1) +
#  geom_text(data = slope_labels_biomass,
 #           aes(x = x, y = y, label = label),
 #           hjust = 0, size = 4, color = "black", check_overlap = TRUE) +
  facet_wrap(~ week_number, nrow = 1, drop = TRUE) +
  guides(alpha = "none", shape = "none", fill = "none") +
  scale_color_brewer(palette="RdBu", direction = -1) +
  scale_shape_manual(values = c("none" = 16, "low" = 17, "high" = 15)) +
  scale_linetype_manual(values = c("none" = "dotdash", "low" = "dashed", "high" = "solid")) +
  labs(
    x = expression(ln~Autotroph~biomass~"(µg C L"^"-1"~")"),
    y = expression(ln~GPP~"(µmol L"^"-1"~"h"^"-1"*")"),
    color = "Wattage", linetype = "Dispersal"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position   = "bottom",
    strip.background  = element_rect(fill = "white", color = "black"),
    panel.background  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA),
    plot.background   = element_blank(),
    panel.grid        = element_blank(),
    strip.text        = element_text(size = 16, family = "Helvetica"),
    axis.title.x      = element_text(size = 16, family = "Helvetica", margin = margin(t = 11)),
    axis.title.y      = element_text(size = 16, family = "Helvetica", margin = margin(r = 11)),
    axis.text.x       = element_text(size = 15, family = "Helvetica"),
    axis.text.y       = element_text(size = 15, family = "Helvetica"),
    legend.text       = element_text(size = 13, family = "Helvetica"),
    legend.title      = element_text(size = 14, family = "Helvetica"),
    legend.margin     = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width  = unit(1, "cm")
  ) +
  guides(
    color    = guide_legend(order = 1, override.aes = list(size = 3.3)),
    linetype = guide_legend(order = 1, override.aes = list(size = 3.3)),
    shape    = guide_legend(order = 1, override.aes = list(size = 3.3))
  )

biomass_glmm


ggsave(biomass_glmm, file="biom_v_temp.pdf",width=10.95,height = 3.8)

## representing param estimates
sdy <- sd(exp_data$ln_GPP_umol_L_h, na.rm = TRUE)
thr <- 0.20 * sdy   # ROPE = [-thr, +thr]

hdi_in_rope_pct <- function(x, thr, ci = 0.89){
  h <- bayestestR::hdi(x, ci = ci)
  lo <- h$CI_low[1]; hi <- h$CI_high[1]
  overlap <- max(0, min(hi,  thr) - max(lo, -thr))
  len     <- max(hi - lo, 1e-12)
  100 * overlap / len
}

summaries <- map_dfr(names(fits), function(ww) {
  smry <- as.data.frame(summary(fits[[ww]])$fixed)
  smry$parameter <- rownames(smry)
  tibble(
    week = ww,
    parameter = smry$parameter,
    estimate = smry$Estimate,
    lower = smry$`l-95% CI`,
    upper = smry$`u-95% CI`
  )
}) %>%
  mutate(
    week = factor(week,
                  levels = c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12"),
                  ordered = TRUE),
    parameter = factor(parameter,
                       levels = c("Intercept","inv_kTc","ln_fluoro_ug_C_L",
                                  "inv_kTc:dispersallow","inv_kTc:dispersalhigh"),
                       labels = c("Intercept",
                                  "Temperature (1/kTc)",
                                  "ln Biomass",
                                  "Temp × Dispersal (low)",
                                  "Temp × Dispersal (high)"))
  )

# Compute the % of 89% HDI in ROPE from draws, per week and parameter 
param_keys <- tibble::tibble(
  param_raw = c("b_Intercept",
                "b_inv_kTc",
                "b_ln_fluoro_ug_C_L",
                "b_inv_kTc:dispersallow",
                "b_inv_kTc:dispersalhigh"),
  parameter = factor(c("Intercept",
                       "Temperature (1/kTc)",
                       "ln Biomass",
                       "Temp × Dispersal (low)",
                       "Temp × Dispersal (high)"),
                     levels = c("Intercept","Temperature (1/kTc)","ln Biomass",
                                "Temp × Dispersal (low)","Temp × Dispersal (high)"))
)

hdi_rope_df <- map_dfr(names(fits), function(ww){
  dr <- posterior::as_draws_df(fits[[ww]])
  keep <- param_keys$param_raw[param_keys$param_raw %in% names(dr)]
  if (length(keep) == 0) return(tibble())
  tmp <- lapply(keep, function(p) {
    vals <- as.numeric(dr[[p]])
    tibble(param_raw = p,
           hdi_in_rope = hdi_in_rope_pct(vals, thr, ci = 0.89))
  }) %>% bind_rows()
  tmp$week <- ww
  tmp
}) %>%
  left_join(param_keys, by = "param_raw") %>%
  mutate(
    week = factor(week,
                  levels = c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12"),
                  ordered = TRUE),
    sig89 = ifelse(hdi_in_rope <= 5, "nonzero", "overlaps")
  ) %>%
  dplyr::select(week, parameter, hdi_in_rope, sig89)

summaries89 <- summaries %>%
  left_join(hdi_rope_df, by = c("week","parameter")) %>%
  mutate(sig89 = tidyr::replace_na(sig89, "overlaps"))

print(summaries89,n=30)

figS6_params <- ggplot(summaries89, aes(x = estimate, y = fct_rev(week), color = sig89)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = sig89), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ parameter, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("nonzero" = "red", "overlaps" = "grey50")) +
  labs(x = "Posterior estimate ±95% CI", y = "Week") +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14, family = "Helvetica"),
    axis.text  = element_text(size = 12, family = "Helvetica"),
    axis.title = element_text(size = 14, family = "Helvetica"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )

supp_fig_params # representing params in supp figure
ggsave(supp_fig_params, file = "gpp_v_temp_mod_params.pdf",width=13,height=3)


