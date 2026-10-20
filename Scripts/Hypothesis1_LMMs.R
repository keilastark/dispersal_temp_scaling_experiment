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

#setwd()
exp_data<-read.csv("metacom_exp_data_final.csv")

exp_data$metacommunity <- factor(exp_data$metacommunity)
exp_data$wattage <- factor(exp_data$wattage, levels=c("0","100","200","300"))
exp_data$dispersal <- factor(exp_data$dispersal, levels=c("none","low","high"))
exp_data$tank <- factor(exp_data$tank)
exp_data$sample_date <- as.Date(exp_data$sample_date, format = "%m/%d/%y") 
exp_data$sample_date <- as.character(exp_data$sample_date)
#unique(exp_data$week_number)

dat2 <- exp_data %>%
  mutate(
    week_number   = factor(
      week_number,
      levels = c("2","4","6","8","10","12"),
      labels = c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12"),
      ordered = TRUE
    ),
    dispersal     = factor(dispersal, levels = c("none","low","high")),
    metacommunity = factor(metacommunity),
    tank          = factor(tank),
    ln_fluoro_ug_C_L = scale(ln_fluoro_ug_C_L) # fluoroprobe-measured chlorophyll, ug C/L
  ) %>%
  mutate(
    dispersal_plot = factor(dispersal, levels = c("none","low","high")),
    week_col_group = if_else(
      week_number %in% c("Week 8"), # strictly for plotting purposes 
      as.character(dispersal_plot),  
      "other"                       
    )
  )

WLVL <- c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12")

priors <- c(
  prior(normal(1.5,5),class = "Intercept"), 
  prior(normal(0.65, 0.28),class = "b", coef = "inv_kTc"), # metabolic theory prior for GPP temperature dependence (E_a)
  prior(normal(0, 0.5), class = "b", coef = "inv_kTc:dispersallow"), # uninformative prior for dispersal effect on E_GPP
  prior(normal(0, 0.5), class = "b", coef = "inv_kTc:dispersalhigh"), 
  prior(exponential(1),   class = "sigma"),
  prior(exponential(1),   class = "sd")
)

# function to fit LMM to each week separately
fit_one_week <- function(d) {
  dat_d <- dat2 %>%
    filter(week_number == d) %>%
    tidyr::drop_na(ln_gpp_corr, inv_kTc,  dispersal, metacommunity, tank) # ln_gpp_corr = ln mass-corrected GPP
  
  stopifnot(nrow(dat_d) >= 8)
  
  brm(
    ln_gpp_corr ~ inv_kTc+inv_kTc:dispersal + (1 |metacommunity), 
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

weeks <- levels(dat2$week_number)
fits  <- setNames(map(weeks, fit_one_week), weeks)
# print fits to view coeff estimates and 95% credible intervals for each week-specific fit


#  ROPE ANALYSIS ON PARAMETER SIGNIFICANCE --------------------------------

sdy <- sd(exp_data$ln_gpp_corr, na.rm = TRUE) # estimate standard deviation in response variable
thr <- 0.20 * sdy   # Set 20% Region of Practical Equivalence (ROPE)

# Function for calculating 89% credible interval on each posterior
hdi_in_rope_pct <- function(x, thr, ci = 0.89){
  h <- bayestestR::hdi(x, ci = ci)
  lo <- h$CI_low[1]; hi <- h$CI_high[1]
  overlap <- max(0, min(hi,  thr) - max(lo, -thr))
  len     <- max(hi - lo, 1e-12)
  100 * overlap / len
}

# Estimate for all parameters and summarize in nicely labelled object
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
                       levels = c("Intercept","inv_kTc",
                                  "inv_kTc:dispersallow","inv_kTc:dispersalhigh"),
                       labels = c("Intercept",
                                  "Temperature (1/kTc)",
                                  "Temp × Dispersal (low)",
                                  "Temp × Dispersal (high)"))
  )

# Compute the % of 89% HDI in ROPE from draws, per week and parameter 
param_keys <- tibble::tibble(
  param_raw = c("b_Intercept",
                "b_inv_kTc",
                "b_inv_kTc:dispersallow",
                "b_inv_kTc:dispersalhigh"),
  parameter = factor(c("Intercept",
                       "Temperature (1/kTc)",
                       "Temp × Dispersal (low)",
                       "Temp × Dispersal (high)"),
                     levels = c("Intercept","Temperature (1/kTc)",
                                "Temp × Dispersal (low)","Temp × Dispersal (high)"))
)

# summarize whether each posterior has nonzero effect based on ROPE criterion
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

hdi_bounds_df <- map_dfr(names(fits), function(ww){
  dr <- posterior::as_draws_df(fits[[ww]])
  keep <- param_keys$param_raw[param_keys$param_raw %in% names(dr)]
  if (length(keep) == 0) return(tibble())
  
  tmp <- lapply(keep, function(p) {
    vals <- as.numeric(dr[[p]])
    hd   <- bayestestR::hdi(vals, ci = 0.89)
    tibble(
      param_raw = p,
      lower = hd$CI_low[1],
      upper = hd$CI_high[1]
    )
  }) %>% bind_rows()
  
  tmp$week <- ww
  tmp
}) %>%
  left_join(param_keys, by = "param_raw") %>%
  mutate(
    week = factor(week,
                  levels  = c("Week 2","Week 4","Week 6","Week 8","Week 10","Week 12"),
                  ordered = TRUE)
  ) %>%
  dplyr::select(week, parameter, lower, upper)

summaries89 <- summaries %>%
  dplyr::select(-lower, -upper) %>%           # drop 95% CI columns, not using
  left_join(hdi_bounds_df,  by = c("week","parameter")) %>%  # add 89% HDI
  left_join(hdi_rope_df,    by = c("week","parameter")) %>%  # add ROPE info
  mutate(sig89 = tidyr::replace_na(sig89, "overlaps"))

print(summaries89, n = 36)


###### Fig 3 Plot Code ######## ----------------------

# extract the temperature dependence slopes from each dispersal level 
.temp_slopes_from_fit <- function(fit) {
  d   <- posterior::as_draws_df(fit)
  bT  <- d$`b_inv_kTc`
  bTL <- if ("b_inv_kTc:dispersallow"  %in% names(d)) d$`b_inv_kTc:dispersallow`  else 0
  bTH <- if ("b_inv_kTc:dispersalhigh" %in% names(d)) d$`b_inv_kTc:dispersalhigh` else 0
  
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

# function to plot predicted regression lines from GLMM fits
pred_lines <- map_dfr(names(fits), function(ww) {
  fit  <- fits[[ww]]
  dat_w <- dat2 %>% filter(week_number == ww) %>%
    drop_na(inv_kTc,dispersal)
  
  xr_w <- range(dat_w$inv_kTc, na.rm = TRUE)
  inv_seq <- seq(xr_w[1] - 0.2, xr_w[2] + 0.2, length.out = 160)  
  
  newd <- expand.grid(inv_kTc = inv_seq, dispersal = levels(dat_w$dispersal)) %>%
    tibble::as_tibble() %>%
    mutate(
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

pred_lines_plot <- pred_lines %>%
  mutate(
    week_col_group = if_else(
      week_number %in% c("Week 8"), # manually color week 8 based on GLMM result
      dispersal_plot,     
      "other"
    )
  )

disp_levels <- c("none","low","high")

# positioning dispersal treatment-specific slope labels
xr <- range(dat2$inv_kTc, na.rm = TRUE)
yr <- c(-1, 0.6)

x_pad <- 0.006 * diff(xr)
y_pad <- 0.08 * diff(yr)
y_gap <- 0.13 * diff(yr)

lab_pos <- tibble(
  dispersal = factor(disp_levels, levels = disp_levels),
  x = xr[1] + x_pad - 0.15,
  y = (yr[2] - y_pad) + 0.3 - (seq_along(disp_levels) - 1) * y_gap - 0.3
)

# Print dispersal treatment-specific slopes from fits onto each facet:
slopes_lab <- slopes_by_week %>%
  left_join(lab_pos, by = "dispersal") %>%
  mutate(
    lab_str = sprintf("%s: %.2f eV", dispersal, slope_mean),
    is_bold = 
      (week_number == "Week 8"  & dispersal == "high"),
    label   = if_else(
      is_bold,
      sprintf("bold('%s')", lab_str),
      sprintf("'%s'", lab_str)
    )
  )

### Plot Manuscript Figure 3(a)
fig3 <- ggplot() +
  geom_point(
    data = dat2,
    aes(inv_kTc, ln_gpp_corr,
        color  = week_col_group,
        shape  = dispersal_plot),
    size = 2.3, alpha = 0.4
  ) +
  geom_ribbon(
    data = pred_lines_plot,
    aes(inv_kTc, ymin = y_lo, ymax = y_hi, fill = week_col_group),
    alpha = 0.07
  ) +
  geom_line(
    data = pred_lines_plot,
    aes(inv_kTc, y = y_mean,
        color    = week_col_group,
        linetype = dispersal_plot),
    linewidth = 1.4
  ) +
  geom_text(
    data = slopes_lab,
    aes(x = x, y = y, label = label),
    hjust = 0, size = 4.5, color = "black",
    parse = TRUE
  ) +
  facet_wrap(~ week_number, nrow = 1, drop = TRUE) +
  guides(alpha = "none", fill = "none") +
  scale_color_manual(
    name   = "Dispersal",
    values = c(
      "none"  = "#003f5c",  # colored dispersal levels
      "low"   = "#ff6361",
      "high"  = "#ffa600",
      "other" = "grey40"    # all non–Week 8 & 10 points/lines
    ),
    breaks = c("none","low","high"),
    labels = c("None","Low","High")
  ) +
  scale_fill_manual(
    values = c(
      "none"  = "#003f5c",
      "low"   = "#ff6361",
      "high"  = "#ffa600",
      "other" = "grey40"
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    name   = "Dispersal",
    values = c("none"=19,"low"=17,"high"=15),
    labels = c("none"="None","low"="Low","high"="High")
  )+
  scale_linetype_manual(
    name   = "Dispersal",
    values = c("none"="dotdash","low"="dashed","high"="solid"),
    labels = c("none"="None","low"="Low","high"="High")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = c("dotdash","dashed","solid"),
        fill     = c("#003f5c", "#ff6361", "#ffa600")
      ),
      order = 1
    ),
    linetype = "none",
    fill     = "none"
  ) +
  labs(
    x = expression("Standardized temperature, " * 1 / italic(k)*italic(T)[c] - 1 / italic(k)*italic(T)*" (1/eV)"),
    y = expression(ln~GPP~"(µmol g C"^"-1"~"h"^"-1"*")"),
    color = "Dispersal", linetype = "Dispersal", shape = "Dispersal"
  ) +scale_x_continuous(breaks = c(-1, 0, 1),
                        labels = c("-1", "0", "1"))+
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
    axis.text.x      = element_text(size = 13, family = "Helvetica"),
    axis.text.y      = element_text(size = 13, family = "Helvetica"),
    legend.text      = element_text(size = 14, family = "Helvetica"),
    legend.title     = element_text(size = 14, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(1.5, "cm")
  ) +ylim(-1,0.5)


print(fig3) 
ggsave(fig3, file="fig3_gpp_v_temp.pdf",width=14.6,height = 3.7)


###

col_nonzero <- "#C62828"
col_other   <- "grey50"

# shared theme for both figures
forest_theme <- theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text       = element_text(size = 11, family = "Helvetica"),
    panel.border     = element_rect(colour = "black", fill = NA),
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.title.x     = element_text(size = 12, family = "Helvetica",
                                    margin = margin(t = 8)),
    axis.title.y     = element_text(size = 12, family = "Helvetica",
                                    margin = margin(r = 8)),
    axis.text        = element_text(size = 10, family = "Helvetica"),
    plot.background  = element_blank(),
    legend.position  = "none"
  )

perweek_rows <- purrr::map_dfr(names(fits), function(ww) {
  d    <- posterior::as_draws_df(fits[[ww]])
  get_d <- function(v) if (v %in% names(d)) as.numeric(d[[v]]) else rep(0, nrow(d))
  
  int_d   <- get_d("b_Intercept")
  temp_d  <- get_d("b_inv_kTc")
  tlow_d  <- get_d("b_inv_kTc:dispersallow")
  thigh_d <- get_d("b_inv_kTc:dispersalhigh")
  
  make_row <- function(draws, param) {
    tibble(
      Week      = factor(ww, levels = rev(WLVL)),
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
  
  bind_rows(
    make_row(int_d,   "Intercept"),
    make_row(temp_d,  "Temperature (1/kTc)"),
    make_row(tlow_d,  "Temp \u00d7 Dispersal (low)"),
    make_row(thigh_d, "Temp \u00d7 Dispersal (high)")
  )
})

panel_rename_pw4 <- c(
  "Temperature (1/kTc)"          = "Temp",
  "Temp \u00d7 Dispersal (low)"  = "Temp~'\u00d7'~Disp~(L~vs~N)",
  "Temp \u00d7 Dispersal (high)" = "Temp~'\u00d7'~Disp~(H~vs~N)"
)

panels_keep4 <- c(
  "Intercept",
  "Temp",
  "Temp~'\u00d7'~Disp~(L~vs~N)",
  "Temp~'\u00d7'~Disp~(H~vs~N)"
)

perweek_sub4 <- perweek_rows %>%
  mutate(Parameter = dplyr::recode(as.character(Parameter), !!!panel_rename_pw4)) %>%
  filter(Parameter %in% panels_keep4) %>%
  mutate(
    Parameter  = factor(Parameter, levels = panels_keep4),
    color_flag = if_else(Decision == "nonzero", "nonzero", "other")
  )

rope_rects4 <- tibble(
  Parameter = factor(panels_keep4, levels = panels_keep4),
  xmin      = -thr,
  xmax      =  thr
)

perweek_sub4 <- perweek_rows %>%
  mutate(Parameter = dplyr::recode(as.character(Parameter), !!!panel_rename_pw4)) %>%
  filter(Parameter %in% panels_keep4) %>%
  mutate(
    Parameter  = factor(Parameter, levels = panels_keep4),
    color_flag = if_else(Decision == "nonzero", "nonzero", "other")
  )

rope_rects4 <- tibble(
  Parameter = factor(panels_keep4, levels = panels_keep4),
  xmin      = -thr,
  xmax      =  thr
)

fig3_forest <- ggplot(
  perweek_sub4,
  aes(y = Week, x = Mean, colour = color_flag)
) +
  geom_rect(
    data        = rope_rects4,
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
    height = 0.25, linewidth = 0.7
  ) +
  geom_point(size = 2.8) +
  scale_colour_manual(
    values = c("nonzero" = col_nonzero, "other" = col_other)
  ) +
  facet_wrap(
    ~Parameter,
    nrow     = 1,
    scales   = "free_x",
    labeller = label_parsed   # <-- this renders the plotmath strings
  ) +  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(
    y = "",
    x = "Posterior mean with 89% HDI"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 13, family = "Helvetica"),
    axis.title.x     = element_text(size = 16, family = "Helvetica",
                                    margin = margin(t = 11)),
    axis.title.y     = element_text(size = 16, family = "Helvetica",
                                    margin = margin(r = 11)),
    axis.text.x      = element_text(size = 13, family = "Helvetica"),
    axis.text.y      = element_text(size = 13, family = "Helvetica"),
    panel.spacing = unit(0.8, "cm")
  )

# Build a dummy plot just to extract a clean legend from
legend_data <- tibble(
  label = factor(
    c("Non-zero effect", "Ambiguous effect", "ROPE"),
    levels = c("Non-zero effect", "Ambiguous effect", "ROPE")
  ),
  x = c(1, 1, 1),
  y = c(1, 2, 3)
)

legend_plot <- ggplot(legend_data) +
  # points for nonzero and ambiguous
  geom_point(
    data = legend_data %>% filter(label != "ROPE"),
    aes(x, y, colour = label),
    size = 3
  ) +
  # rectangle for ROPE
  geom_rect(
    data = legend_data %>% filter(label == "ROPE"),
    aes(xmin = 0.8, xmax = 1.2, ymin = 2.7, ymax = 3.3, fill = label),
    colour = NA
  ) +
  scale_colour_manual(
    name   = NULL,
    values = c(
      "Non-zero effect"  = col_nonzero,
      "Ambiguous effect" = col_other
    ),
    guide = guide_legend(order = 1)
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c("ROPE" = "grey85"),
    guide  = guide_legend(order = 2)
  ) +
  theme_void(base_size = 13) +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 12, family = "Helvetica"),
    legend.key       = element_rect(fill = "white", colour = NA),
    legend.spacing.y = unit(0.6, "cm")
  )

forest_legend <- cowplot::get_legend(legend_plot)

# Attach legend to right of forest plot
fig3_forest_with_legend <- cowplot::plot_grid(
  fig3_forest,
  forest_legend,
  nrow       = 1,
  rel_widths = c(1, 0.15)
)

fig3_forest_with_legend
#ggsave(fig3_forest_with_legend, file   = "fig3_forest.pdf", width  = 12.7,height = 2.6)


# Stack

fig_combined <- fig3 / fig3_forest_with_legend +
  plot_layout(heights = c(2, 2.5 )) +
  plot_annotation(tag_levels = "a")
fig_combined # essentially replicates fig3 in ms, however the individual panels were compiled in graphics program for tidiness




