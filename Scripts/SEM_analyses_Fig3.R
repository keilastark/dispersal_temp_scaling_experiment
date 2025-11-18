library(dplyr)
library(tidyr)
library(brms)
library(posterior)
library(bayestestR)
library(purrr)
library(tibble)
library(ggplot2)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# install_cmdstan() #If you want to use cmdstan backend 

exp_data<-read.csv("metacom_exp_data_final.csv")

exp_data$metacommunity <- factor(exp_data$metacommunity)
exp_data$wattage <- factor(exp_data$wattage, levels=c("0","100","200","300"))
exp_data$dispersal <- factor(exp_data$dispersal, levels=c("none","low","high"))
exp_data$tank <- factor(exp_data$tank)
exp_data$sample_date <- as.character(exp_data$sample_date)

df_sem_subset <- exp_data %>%
  dplyr::select(
    ln_GPP_umol_L_h,
    size_PCoA1, comp_PCoA1, size_NMDS1,comp_NMDS1,
    ln_fluoro_ug_C_L, inv_kTc,
    totalN,
    dispersal, metacommunity, tank, ln_zoop_ug_L, week_number,size_NMDS_mean
  ) %>%
  tidyr::drop_na()

df_sem_subset <- df_sem_subset %>%
  mutate(
    inv_kTc             = inv_kTc,
    size_NMDS        = size_NMDS1,
    comp_NMDS        = comp_NMDS1,
    ln_fluoro_ug_C_L  = ln_fluoro_ug_C_L,
    ln_totalN                 = log(totalN),
    ln_GPP_umol_L_h   = as.numeric(ln_GPP_umol_L_h),
    dispersal           = factor(dispersal),
    metacommunity       = factor(metacommunity),
    tank                = factor(tank),
    zoop_ug_C_L       = ln_zoop_ug_L,
    week_number         = factor(week_number, levels = c("4","8","12"), ordered = FALSE)
  )


# Model formulas & priors (shared across weeks)

bf_size <- bf(size_NMDS ~ dispersal + inv_kTc + zoop_ug_C_L + (1|metacommunity), family = gaussian())
bf_comp <- bf(comp_NMDS ~ dispersal + inv_kTc + zoop_ug_C_L + (1|metacommunity), family = gaussian())
bf_biomass <- bf(ln_fluoro_ug_C_L ~ dispersal + inv_kTc + ln_totalN + zoop_ug_C_L + (1 | metacommunity))
bf_zoop <- bf(zoop_ug_C_L ~ inv_kTc + dispersal + (1|metacommunity), family = gaussian())

bf_gpp_all            <- bf(ln_GPP_umol_L_h ~ size_NMDS + comp_NMDS + ln_fluoro_ug_C_L + inv_kTc + ln_totalN + (1|metacommunity), family = gaussian())


priors_all <- c(
  # Intercepts
  set_prior("normal(0,1)", class = "Intercept", resp = "sizeNMDS"),   
  set_prior("normal(0,1)", class = "Intercept", resp = "compNMDS"),
  set_prior("normal(0,1)", class = "Intercept", resp = "lnfluorougCL"),
  set_prior("normal(0,1)", class = "Intercept", resp = "lnGPPumolLh"),
  set_prior("normal(0,1)", class = "Intercept", resp = "zoopugCL"),
  
  # Size structure (phyto)
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc",       resp = "sizeNMDS"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow",  resp = "sizeNMDS"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh", resp = "sizeNMDS"),
  set_prior("normal(0.5,0.5)", class = "b", coef = "zoop_ug_C_L", resp = "sizeNMDS"),
  
  # Composition (phyto)
  set_prior("normal(0,0.5)", class = "b", coef = "inv_kTc",       resp = "compNMDS"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersallow",  resp = "compNMDS"),
  set_prior("normal(0,0.5)", class = "b", coef = "dispersalhigh", resp = "compNMDS"),
  set_prior("normal(0,0.5)", class = "b", coef = "zoop_ug_C_L",   resp = "compNMDS"),
  
  # Biomass (phyto)
  set_prior("normal(-0.5,0.5)", class = "b", coef = "inv_kTc",       resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)",    class = "b", coef = "dispersallow",  resp = "lnfluorougCL"),
  set_prior("normal(0,0.5)",    class = "b", coef = "dispersalhigh", resp = "lnfluorougCL"),
  set_prior("normal(1,0.5)",    class = "b", coef = "ln_totalN",     resp = "lnfluorougCL"),
  set_prior("normal(-0.5,0.5)", class = "b", coef = "zoop_ug_C_L",   resp = "lnfluorougCL"),
  
  # Zooplankton (biomass)
  set_prior("normal(0.5,0.5)", class = "b", coef = "inv_kTc",       resp = "zoopugCL"),
  set_prior("normal(0,0.5)",   class = "b", coef = "dispersallow",  resp = "zoopugCL"),
  set_prior("normal(0,0.5)",   class = "b", coef = "dispersalhigh", resp = "zoopugCL"),
  
  # GPP 
  set_prior("normal(0.65,0.28)", class = "b", coef = "inv_kTc",       resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)",     class = "b", coef = "ln_totalN",     resp = "lnGPPumolLh"),
  set_prior("normal(0,0.5)",     class = "b", coef = "size_NMDS",     resp = "lnGPPumolLh"), 
  set_prior("normal(0,0.5)",     class = "b", coef = "comp_NMDS",     resp = "lnGPPumolLh"), 
  set_prior("normal(0.75,0.5)",  class = "b", coef = "ln_fluoro_ug_C_L", resp = "lnGPPumolLh"),
  
  # Residual SDs 
  set_prior("exponential(1)", class = "sigma", resp = "sizeNMDS"),
  set_prior("exponential(1)", class = "sigma", resp = "compNMDS"),
  set_prior("exponential(1)", class = "sigma", resp = "lnfluorougCL"),
  set_prior("exponential(1)", class = "sigma", resp = "lnGPPumolLh"),
  set_prior("exponential(1)", class = "sigma", resp = "zoopugCL")
)

## % HDI in ROPE helper function
hdi_in_rope_pct <- function(x, thr, ci = 0.89){
  h <- bayestestR::hdi(x, ci = ci)
  lo <- h$CI_low[1]; hi <- h$CI_high[1]
  overlap <- max(0, min(hi,  thr) - max(lo, -thr))
  len     <- max(hi - lo, 1e-12)
  100 * overlap / len
}

## =========================================================
## WEEK 4
## =========================================================
df_sem_subset_week4 <- df_sem_subset %>% filter(week_number == "4")

fit_all_week4 <- brm(
  bf_size + bf_comp + bf_biomass + bf_gpp_all + bf_zoop + set_rescor(FALSE),
  data = df_sem_subset_week4, prior = priors_all,
  backend = "cmdstanr", chains = 4, cores = 4, iter = 5000,
   save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15), seed = 2025
)

draws_week4 <- posterior::as_draws_df(fit_all_week4)
get_coef_week4 <- function(var) as.numeric(draws_week4[[var]])

# Direct paths to GPP
size_gpp_week4   <- get_coef_week4("b_lnGPPumolLh_size_NMDS")
comp_gpp_week4   <- get_coef_week4("b_lnGPPumolLh_comp_NMDS")
phyto_gpp_week4  <- get_coef_week4("b_lnGPPumolLh_ln_fluoro_ug_C_L")
temp_gpp_week4   <- get_coef_week4("b_lnGPPumolLh_inv_kTc")
N_gpp_week4      <- get_coef_week4("b_lnGPPumolLh_ln_totalN")

# Temp paths into mediators 
temp_size_week4  <- get_coef_week4("b_size_NMDS_inv_kTc")
temp_comp_week4  <- get_coef_week4("b_compNMDS_inv_kTc")
temp_phyto_week4 <- get_coef_week4("b_lnfluorougCL_inv_kTc")
temp_zoop_week4  <- get_coef_week4("b_zoopugCL_inv_kTc")

#Nitrogen direct effect
N_phyto_week4    <- get_coef_week4("b_lnfluorougCL_ln_totalN")

# Dispersal paths into mediators
low_size_week4   <- get_coef_week4("b_size_NMDS_dispersallow")
high_size_week4  <- get_coef_week4("b_size_NMDS_dispersalhigh")
low_comp_week4   <- get_coef_week4("b_compNMDS_dispersallow")
high_comp_week4  <- get_coef_week4("b_compNMDS_dispersalhigh")
low_phyto_week4  <- get_coef_week4("b_lnfluorougCL_dispersallow")
high_phyto_week4 <- get_coef_week4("b_lnfluorougCL_dispersalhigh")

low_zoop_week4   <- get_coef_week4("b_zoopugCL_dispersallow")
high_zoop_week4  <- get_coef_week4("b_zoopugCL_dispersalhigh")

zoop_size_week4  <- get_coef_week4("b_size_NMDS_zoop_ug_C_L")
zoop_comp_week4  <- get_coef_week4("b_compNMDS_zoop_ug_C_L")
zoop_phyto_week4 <- get_coef_week4("b_lnfluorougCL_zoop_ug_C_L")

# Totals / indirects
ind_temp_size_week4  <- temp_size_week4 * size_gpp_week4
ind_temp_comp_week4  <- temp_comp_week4 * comp_gpp_week4
ind_temp_phyto_week4 <- temp_phyto_week4 * phyto_gpp_week4
temp_total_week4     <- temp_gpp_week4 + ind_temp_size_week4 + ind_temp_comp_week4 + ind_temp_phyto_week4

ind_N_phyto_week4    <- N_phyto_week4 * phyto_gpp_week4
N_total_week4        <- N_gpp_week4 + ind_N_phyto_week4

# dispersal indirects via (size, comp, phyto)
ind_low_week4  <- low_size_week4  * size_gpp_week4 + low_comp_week4  * comp_gpp_week4 + low_phyto_week4  * phyto_gpp_week4
ind_high_week4 <- high_size_week4 * size_gpp_week4 + high_comp_week4 * comp_gpp_week4 + high_phyto_week4 * phyto_gpp_week4

# Table of edge posteriors
edges_draws_week4 <- tibble::tribble(
  ~from,            ~to,              ~vec,
  "Temp",           "GPP",            temp_gpp_week4,
  "Total N",        "GPP",            N_gpp_week4,
  "Size structure", "GPP",            size_gpp_week4,
  "Composition",    "GPP",            comp_gpp_week4,
  "Phyto biomass",  "GPP",            phyto_gpp_week4,
  "Temp",           "Zoop biomass",   temp_zoop_week4,
  "Temp",           "Size structure", temp_size_week4,
  "Temp",           "Composition",    temp_comp_week4,
  "Temp",           "Phyto biomass",  temp_phyto_week4,
  "Total N",        "Phyto biomass",  N_phyto_week4,
  "Dispersal: low", "Zoop biomass",   low_zoop_week4,
  "Dispersal: high","Zoop biomass",   high_zoop_week4,
  "Dispersal: low", "Size structure", low_size_week4,
  "Dispersal: high","Size structure", high_size_week4,
  "Dispersal: low", "Composition",    low_comp_week4,
  "Dispersal: high","Composition",    high_comp_week4,
  "Dispersal: low", "Phyto biomass",  low_phyto_week4,
  "Dispersal: high","Phyto biomass",  high_phyto_week4,
  "Zoop biomass",   "Size structure", zoop_size_week4,
  "Zoop biomass",   "Composition",    zoop_comp_week4,
  "Zoop biomass",   "Phyto biomass",  zoop_phyto_week4
)


sdy_week4 <- sd(df_sem_subset_week4$ln_GPP_umol_L_h, na.rm = TRUE)
thr_week4 <- 0.20*sdy_week4

# flagging whether HDI in ROPE
edge_summ_week4 <- edges_draws_week4 %>%
  mutate(
    Mean      = purrr::map_dbl(vec, mean),
    L89       = purrr::map_dbl(vec, ~quantile(.x, 0.055)),
    U89       = purrr::map_dbl(vec, ~quantile(.x, 0.945)),
    HDIinROPE = purrr::map_dbl(vec, ~hdi_in_rope_pct(.x, thr_week4,  ci = 0.89)),
    decision  = case_when(HDIinROPE >= 95 ~ "null",
                          HDIinROPE <=  5 ~ "nonzero",
                          TRUE            ~ "ambiguous")
  )

# Following section of code setting aesthetics for DAG plot:
style_df_week4  <- edge_summ_week4 %>%  #
  transmute(from, to, Mean, L89, U89, HDIinROPE, decision,
            sign = ifelse(Mean >= 0, "pos", "neg"),
            is_sig = decision == "nonzero",
            color = case_when(is_sig & sign=="pos" ~ "#222222",
                              is_sig & sign=="neg" ~ "#C62828",
                              TRUE                 ~ "#878787"),
            style = ifelse(is_sig, "solid", "dotted"),
            penwidth = ifelse(is_sig, 5, 4),
            arrowsize = ifelse(is_sig, 1.2, 0.9))
tg4 <- dplyr::filter(style_df_week4, from == "Temp", to == "GPP") %>% dplyr::slice(1)
if (nrow(tg4) != 1) stop("Temp->GPP style not found (week 4)")

top_nodes <- c("Temp", "Dispersal: low", "Dispersal: high")
mid2_nodes <- c("Zoop biomass", "Total N")
mid3_nodes <- c("Phyto biomass", "Composition", "Size structure")
bot_nodes <- c("GPP")

dot_header <- paste(
  'digraph sem {',
  '  graph [layout=dot, rankdir=TB, newrank=true,',
  '         splines=curved, overlap=false,',
  '         nodesep="0.8", ranksep="1.8",',
  '         margin=0.05, pad="0.1,0.1", center=true];',
  '  node  [shape=ellipse, style="filled", fillcolor="white", color="black",',
  '         fontname="Helvetica", fontsize=20, fixedsize=true, width=2, height=1];',
  '  edge  [fontname="Helvetica", fontsize=14, color="#888888", arrowsize=0.9, labelfloat=false];',
  sep = "\n"
)

all_nodes_week4 <- unique(c(style_df_week4$from, style_df_week4$to, top_nodes, mid2_nodes, mid3_nodes, bot_nodes))
node_lines_week4 <- vapply(all_nodes_week4, function(lbl) sprintf('  "%s" [shape=ellipse];', lbl), character(1))

rank_top <- paste(
  '  { rank=same;',
  '    Ctop [shape=point, width=0, label="", style=invis];',
  '    SpacerTopL [shape=point, label="", style=invis, width=0.6, height=0.01];',
  '    "Temp";',
  '    "Dispersal: low";',
  '    "Dispersal: high";',
  '  }', sep = "\n"
)
rank_mid2_redo <- paste(
  '  { rank=same;',
  '    Cmid2 [shape=point, width=0, label="", style=invis];',
  '    "Zoop biomass";',
  '    SpacerN [shape=point, label="", style=invis, width=0.6, height=0.01];',
  '    "Total N";',
  '  }', sep = "\n"
)
rank_mid3 <- paste(
  '  { rank=same;',
  '    Cmid3 [shape=point, width=0, label="", style=invis];',
  '    "Size structure";',
  '    Spacer3a [shape=point, label="", style=invis, width=0.6, height=0.01];',
  '    "Composition";',
  '    Spacer3b [shape=point, label="", style=invis, width=0.6, height=0.01];',
  '    "Phyto biomass";',
  '  }', sep = "\n"
)
rank_bot <- paste(
  '  { rank=same;',
  '    Cbot [shape=point, width=0, label="", style=invis];',
  '    "GPP";',
  '  }', sep = "\n"
)
horiz_top <- paste(
  '  { rank=same;',
  '    SpacerTopL -> "Temp" -> "Dispersal: low" -> "Dispersal: high" [style=invis, weight=80];',
  '  }', sep = "\n"
)
horiz_mid2 <- paste(
  '  { rank=same;',
  '    "Zoop biomass" -> SpacerN -> "Total N" [style=invis, weight=70];',
  '  }', sep = "\n"
)
horiz_mid3 <- paste(
  '  { rank=same;',
  '    "Size structure" -> Spacer3a -> "Composition" -> Spacer3b -> "Phyto biomass" [style=invis, weight=85];',
  '  }', sep = "\n"
)
keep_N_right <- paste(
  '  "Dispersal: high" -> "Total N" [style=invis, weight=51];',
  '  "Phyto biomass"   -> "Total N" [style=invis, weight=47];', sep = "\n"
)
center_helps <- paste(
  '  "Dispersal: low"  -> "Zoop biomass" [style=invis, weight=20];',
  '  "Dispersal: high" -> "Zoop biomass" [style=invis, weight=20];',
  '  "Size structure"  -> "Zoop biomass" [style=invis, weight=20];',
  '  "Phyto biomass"   -> "Zoop biomass" [style=invis, weight=20];', sep = "\n"
)
center_lines <- c(
  '  Ctop -> Cmid2 -> Cmid3 -> Cbot [style=invis, weight=200];',
  paste(sprintf('  Ctop  -> "%s" [style=invis, constraint=false, weight=80];', top_nodes), collapse = "\n"),
  paste(sprintf('  Cmid2 -> "%s" [style=invis, constraint=false, weight=80];', mid2_nodes), collapse = "\n"),
  paste(sprintf('  Cmid3 -> "%s" [style=invis, constraint=false, weight=80];', mid3_nodes), collapse = "\n"),
  '  Cbot  -> "GPP" [style=invis, constraint=false, weight=80];'
)

style_df_no_TG_week4 <- subset(style_df_week4, !(from == "Temp" & to == "GPP"))
edge_lines_week4 <- vapply(
  seq_len(nrow(style_df_no_TG_week4)),
  function(i) sprintf(
    '  "%s" -> "%s" [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s];',
    style_df_no_TG_week4$from[i], style_df_no_TG_week4$to[i],
    style_df_no_TG_week4$color[i], style_df_no_TG_week4$penwidth[i],
    style_df_no_TG_week4$arrowsize[i], style_df_no_TG_week4$style[i]
  ), character(1)
)
temp_gpp_left_week4 <- sprintf(
  '  "Temp":s -> "GPP":nw [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s, weight=70, minlen=2, constraint=false];',
  tg4$color, tg4$penwidth, tg4$arrowsize, tg4$style
)

newdot_week4 <- paste(
  dot_header,
  paste(node_lines_week4, collapse = "\n"),
  rank_top, rank_mid2_redo, rank_mid3, rank_bot,
  horiz_top, horiz_mid2, horiz_mid3, keep_N_right, center_helps,
  paste(edge_lines_week4, collapse = "\n"),
  temp_gpp_left_week4,
  paste(center_lines, collapse = "\n"),
  "}", sep = "\n"
)

SEM_fig_week4 <- DiagrammeR::grViz(newdot_week4, width = "100%")
SEM_fig_week4
svg4 <- export_svg(SEM_fig_week4)
# rsvg_pdf(charToRaw(svg4), file = "sem_diagram_week4.pdf", width = 1210, height = 590) 

#ROPE density plot inputs (direct + indirect sets) 
effect_draws_week4 <- tibble::tibble(
  Effect = c("Temp→GPP total","N→GPP total",
             "Low dispersal → GPP (ind.)","High dispersal → GPP (ind.)",
             "Size→GPP","Comp→GPP","Biomass→GPP",
             "Temp→GPP (direct)","Temp→GPP (ind.)"),
  Draws  = list(
    temp_total_week4,
    N_total_week4,
    ind_low_week4, ind_high_week4,
    size_gpp_week4, comp_gpp_week4, phyto_gpp_week4,
    temp_gpp_week4, (ind_temp_size_week4 + ind_temp_comp_week4 + ind_temp_phyto_week4)
  )
) %>% tidyr::unnest_longer(Draws, indices_include = FALSE, values_to = "value")

rope_band_week4 <- data.frame(xmin = -thr_week4, xmax = thr_week4)
lab_vec_week4 <- setNames(
  paste0(unique(effect_draws_week4$Effect)),
  unique(effect_draws_week4$Effect)
)

### ROPE fig as in Supp Mat
ROPEs_week4 <- ggplot(effect_draws_week4, aes(x = value)) +
  geom_rect(data = rope_band_week4,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.3, fill = "pink3") +
  geom_density(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~ Effect, scales = "free", ncol = 3,
             labeller = labeller(Effect = lab_vec_week4)) +
  labs(x = "Effect (log-GPP scale)", y = "Density",
       title = "Posterior densities with ROPE (±0.20·SD(GPP)) — Week 4") +
  theme_classic(base_size = 12)

## =========================================================
## WEEK 8  (identical workflow)
## =========================================================
df_sem_subset_week8 <- df_sem_subset %>% filter(week_number == "8")

fit_all_week8 <- brm(
  bf_size + bf_comp + bf_biomass + bf_gpp_all + bf_zoop + set_rescor(FALSE),
  data = df_sem_subset_week8, prior = priors_all,
  chains = 4, cores = 4, iter = 5000,
  backend = "cmdstanr", save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15), seed = 2025
)

sdy_week8 <- sd(df_sem_subset_week8$ln_GPP_umol_L_h, na.rm = TRUE)
thr_week8 <- 0.20*sdy_week8

draws_week8 <- posterior::as_draws_df(fit_all_week8)
get_coef_week8 <- function(var) as.numeric(draws_week8[[var]])

size_gpp_week8   <- get_coef_week8("b_lnGPPumolLh_size_NMDS")
comp_gpp_week8   <- get_coef_week8("b_lnGPPumolLh_comp_NMDS")
phyto_gpp_week8  <- get_coef_week8("b_lnGPPumolLh_ln_fluoro_ug_C_L")
temp_gpp_week8   <- get_coef_week8("b_lnGPPumolLh_inv_kTc")
N_gpp_week8      <- get_coef_week8("b_lnGPPumolLh_ln_totalN")

temp_size_week8  <- get_coef_week8("b_size_NMDS_inv_kTc")
temp_comp_week8  <- get_coef_week8("b_compNMDS_inv_kTc")
temp_phyto_week8 <- get_coef_week8("b_lnfluorougCL_inv_kTc")
temp_zoop_week8  <- get_coef_week8("b_zoopugCL_inv_kTc")

N_phyto_week8    <- get_coef_week8("b_lnfluorougCL_ln_totalN")

low_size_week8   <- get_coef_week8("b_size_NMDS_dispersallow")
high_size_week8  <- get_coef_week8("b_size_NMDS_dispersalhigh")
low_comp_week8   <- get_coef_week8("b_compNMDS_dispersallow")
high_comp_week8  <- get_coef_week8("b_compNMDS_dispersalhigh")
low_phyto_week8  <- get_coef_week8("b_lnfluorougCL_dispersallow")
high_phyto_week8 <- get_coef_week8("b_lnfluorougCL_dispersalhigh")

low_zoop_week8   <- get_coef_week8("b_zoopugCL_dispersallow")
high_zoop_week8  <- get_coef_week8("b_zoopugCL_dispersalhigh")

zoop_size_week8  <- get_coef_week8("b_size_NMDS_zoop_ug_C_L")
zoop_comp_week8  <- get_coef_week8("b_compNMDS_zoop_ug_C_L")
zoop_phyto_week8 <- get_coef_week8("b_lnfluorougCL_zoop_ug_C_L")

ind_temp_size_week8  <- temp_size_week8 * size_gpp_week8
ind_temp_comp_week8  <- temp_comp_week8 * comp_gpp_week8
ind_temp_phyto_week8 <- temp_phyto_week8 * phyto_gpp_week8
temp_total_week8     <- temp_gpp_week8 + ind_temp_size_week8 + ind_temp_comp_week8 + ind_temp_phyto_week8

ind_N_phyto_week8    <- N_phyto_week8 * phyto_gpp_week8
N_total_week8        <- N_gpp_week8 + ind_N_phyto_week8

ind_low_week8  <- low_size_week8  * size_gpp_week8 + low_comp_week8  * comp_gpp_week8 + low_phyto_week8  * phyto_gpp_week8
ind_high_week8 <- high_size_week8 * size_gpp_week8 + high_comp_week8 * comp_gpp_week8 + high_phyto_week8 * phyto_gpp_week8

edges_draws_week8_NMDS <- tibble::tribble(
  ~from,            ~to,              ~vec,
  "Temp",           "GPP",            temp_gpp_week8,
  "Total N",        "GPP",            N_gpp_week8,
  "Size structure", "GPP",            size_gpp_week8,
  "Composition",    "GPP",            comp_gpp_week8,
  "Phyto biomass",  "GPP",            phyto_gpp_week8,
  "Temp",           "Zoop biomass",   temp_zoop_week8,
  "Temp",           "Size structure", temp_size_week8,
  "Temp",           "Composition",    temp_comp_week8,
  "Temp",           "Phyto biomass",  temp_phyto_week8,
  "Total N",        "Phyto biomass",  N_phyto_week8,
  "Dispersal: low", "Zoop biomass",   low_zoop_week8,
  "Dispersal: high","Zoop biomass",   high_zoop_week8,
  "Dispersal: low", "Size structure", low_size_week8,
  "Dispersal: high","Size structure", high_size_week8,
  "Dispersal: low", "Composition",    low_comp_week8,
  "Dispersal: high","Composition",    high_comp_week8,
  "Dispersal: low", "Phyto biomass",  low_phyto_week8,
  "Dispersal: high","Phyto biomass",  high_phyto_week8,
  "Zoop biomass",   "Size structure", zoop_size_week8,
  "Zoop biomass",   "Composition",    zoop_comp_week8,
  "Zoop biomass",   "Phyto biomass",  zoop_phyto_week8
)

edge_summ_week8 <- edges_draws_week8 %>%
      mutate(
        Mean      = purrr::map_dbl(vec, mean),
        L89       = purrr::map_dbl(vec, ~quantile(.x, 0.055)),
        U89       = purrr::map_dbl(vec, ~quantile(.x, 0.945)),
        HDIinROPE = purrr::map_dbl(vec, ~hdi_in_rope_pct(.x, thr_week8,  ci = 0.89)),
        decision  = case_when(HDIinROPE >= 95 ~ "null",
                              HDIinROPE <=  5 ~ "nonzero",
                              TRUE            ~ "ambiguous")
      )
    
style_df_week8  <- edge_summ_week8 %>%  # <- was edge_summ_week8
  transmute(from, to, Mean, L89, U89, HDIinROPE, decision,
            sign = ifelse(Mean >= 0, "pos", "neg"),
            is_sig = decision == "nonzero",
            color = case_when(is_sig & sign=="pos" ~ "#222222",
                              is_sig & sign=="neg" ~ "#C62828",
                              TRUE                 ~ "#878787"),
            style = ifelse(is_sig, "solid", "dotted"),
            penwidth = ifelse(is_sig, 5, 4),
            arrowsize = ifelse(is_sig, 1.2, 0.9))

  
tg8 <- dplyr::filter(style_df_week8, from == "Temp", to == "GPP") %>% dplyr::slice(1)
if (nrow(tg8) != 1) stop("Temp->GPP style not found (week 8)")

all_nodes_week8 <- unique(c(style_df_week8$from, style_df_week8$to, top_nodes, mid2_nodes, mid3_nodes, bot_nodes))
node_lines_week8 <- vapply(all_nodes_week8, function(lbl) sprintf('  "%s" [shape=ellipse];', lbl), character(1))

edge_lines_week8 <- vapply(
  seq_len(nrow(subset(style_df_week8, !(from == "Temp" & to == "GPP")))),
  function(i) {
    sdf <- subset(style_df_week8, !(from == "Temp" & to == "GPP"))
    sprintf('  "%s" -> "%s" [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s];',
            sdf$from[i], sdf$to[i], sdf$color[i], sdf$penwidth[i], sdf$arrowsize[i], sdf$style[i])
  }, character(1)
)
temp_gpp_left_week8 <- sprintf(
  '  "Temp":s -> "GPP":nw [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s, weight=70, minlen=2, constraint=false];',
  tg8$color, tg8$penwidth, tg8$arrowsize, tg8$style
)

newdot_week8 <- paste(
  dot_header,
  paste(node_lines_week8, collapse = "\n"),
  rank_top, rank_mid2_redo, rank_mid3, rank_bot,
  horiz_top, horiz_mid2, horiz_mid3, keep_N_right, center_helps,
  paste(edge_lines_week8, collapse = "\n"),
  temp_gpp_left_week8,
  paste(center_lines, collapse = "\n"),
  "}", sep = "\n"
)

SEM_fig_week8 <- DiagrammeR::grViz(newdot_week8, width = "100%")
SEM_fig_week8
svg8 <- export_svg(SEM_fig_week8)
#rsvg_pdf(charToRaw(svg8), file = "sem_diagram_week8.pdf", width = 1210, height = 590)

effect_draws_week8 <- tibble::tibble(
  Effect = c("Temp→GPP total","N→GPP total",
             "Low dispersal → GPP (ind.)","High dispersal → GPP (ind.)",
             "Size→GPP","Comp→GPP","Biomass→GPP",
             "Temp→GPP (direct)","Temp→GPP (ind.)"),
  Draws  = list(
    temp_total_week8,
    N_total_week8,
    ind_low_week8, ind_high_week8,
    size_gpp_week8, comp_gpp_week8, phyto_gpp_week8,
    temp_gpp_week8, (ind_temp_size_week8 + ind_temp_comp_week8 + ind_temp_phyto_week8)
  )
) %>% tidyr::unnest_longer(Draws, indices_include = FALSE, values_to = "value")

rope_band_week8 <- data.frame(xmin = -thr_week8, xmax = thr_week8)
lab_vec_week8 <- setNames(paste0(unique(effect_draws_week8$Effect)), unique(effect_draws_week8$Effect))

ROPEs_week8 <- ggplot(effect_draws_week8, aes(x = value)) +
  geom_rect(data = rope_band_week8,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.3, fill = "pink3") +
  geom_density(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~ Effect, scales = "free", ncol = 3,
             labeller = labeller(Effect = lab_vec_week8)) +
  labs(x = "Effect (log-GPP scale)", y = "Density",
       title = "Posterior densities with ROPE (±0.20·SD(GPP)) — Week 8") +
  theme_classic(base_size = 12)


## =========================================================
## WEEK 12  (same thing again)
## =========================================================
df_sem_subset_week12 <- df_sem_subset %>% filter(week_number == "12")

fit_all_week12 <- brm(
  bf_size + bf_comp + bf_biomass + bf_gpp_all + bf_zoop + set_rescor(FALSE),
  data = df_sem_subset_week12, prior = priors_all,
  chains = 4, cores = 4, iter = 4000,
  backend = "cmdstanr", save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15), seed = 2025
)
sdy_week12 <- sd(df_sem_subset_week12$ln_GPP_umol_L_h, na.rm = TRUE)
thr_week12 <- 0.20*sdy_week12

draws_week12 <- posterior::as_draws_df(fit_all_week12)
get_coef_week12 <- function(var) as.numeric(draws_week12[[var]])

size_gpp_week12   <- get_coef_week12("b_lnGPPumolLh_size_NMDS")
comp_gpp_week12   <- get_coef_week12("b_lnGPPumolLh_comp_NMDS")
phyto_gpp_week12  <- get_coef_week12("b_lnGPPumolLh_ln_fluoro_ug_C_L")
temp_gpp_week12   <- get_coef_week12("b_lnGPPumolLh_inv_kTc")
N_gpp_week12      <- get_coef_week12("b_lnGPPumolLh_ln_totalN")

temp_size_week12  <- get_coef_week12("b_size_NMDS_inv_kTc")
temp_comp_week12  <- get_coef_week12("b_compNMDS_inv_kTc")
temp_phyto_week12 <- get_coef_week12("b_lnfluorougCL_inv_kTc")
temp_zoop_week12  <- get_coef_week12("b_zoopugCL_inv_kTc")

N_phyto_week12    <- get_coef_week12("b_lnfluorougCL_ln_totalN")

low_size_week12   <- get_coef_week12("b_size_NMDS_dispersallow")
high_size_week12  <- get_coef_week12("b_size_NMDS_dispersalhigh")
low_comp_week12   <- get_coef_week12("b_compNMDS_dispersallow")
high_comp_week12  <- get_coef_week12("b_compNMDS_dispersalhigh")
low_phyto_week12  <- get_coef_week12("b_lnfluorougCL_dispersallow")
high_phyto_week12 <- get_coef_week12("b_lnfluorougCL_dispersalhigh")


low_zoop_week12   <- get_coef_week12("b_zoopugCL_dispersallow")
high_zoop_week12  <- get_coef_week12("b_zoopugCL_dispersalhigh")

zoop_size_week12  <- get_coef_week12("b_size_NMDS_zoop_ug_C_L")
zoop_comp_week12  <- get_coef_week12("b_compNMDS_zoop_ug_C_L")
zoop_phyto_week12 <- get_coef_week12("b_lnfluorougCL_zoop_ug_C_L")

ind_temp_size_week12  <- temp_size_week12 * size_gpp_week12
ind_temp_comp_week12  <- temp_comp_week12 * comp_gpp_week12
ind_temp_phyto_week12 <- temp_phyto_week12 * phyto_gpp_week12
temp_total_week12     <- temp_gpp_week12 + ind_temp_size_week12 + ind_temp_comp_week12 + ind_temp_phyto_week12

ind_N_phyto_week12    <- N_phyto_week12 * phyto_gpp_week12
N_total_week12        <- N_gpp_week12 + ind_N_phyto_week12

ind_low_week12  <- low_size_week12  * size_gpp_week12 + low_comp_week12  * comp_gpp_week12 + low_phyto_week12  * phyto_gpp_week12
ind_high_week12 <- high_size_week12 * size_gpp_week12 + high_comp_week12 * comp_gpp_week12 + high_phyto_week12 * phyto_gpp_week12

edges_draws_week12 <- tibble::tribble(
  ~from,            ~to,              ~vec,
  "Temp",           "GPP",            temp_gpp_week12,
  "Total N",        "GPP",            N_gpp_week12,
  "Size structure", "GPP",            size_gpp_week12,
  "Composition",    "GPP",            comp_gpp_week12,
  "Phyto biomass",  "GPP",            phyto_gpp_week12,
  "Temp",           "Zoop biomass",   temp_zoop_week12,
  "Temp",           "Size structure", temp_size_week12,
  "Temp",           "Composition",    temp_comp_week12,
  "Temp",           "Phyto biomass",  temp_phyto_week12,
  "Total N",        "Phyto biomass",  N_phyto_week12,
  "Dispersal: low", "Zoop biomass",   low_zoop_week12,
  "Dispersal: high","Zoop biomass",   high_zoop_week12,
  "Dispersal: low", "Size structure", low_size_week12,
  "Dispersal: high","Size structure", high_size_week12,
  "Dispersal: low", "Composition",    low_comp_week12,
  "Dispersal: high","Composition",    high_comp_week12,
  "Dispersal: low", "Phyto biomass",  low_phyto_week12,
  "Dispersal: high","Phyto biomass",  high_phyto_week12,
  "Zoop biomass",   "Size structure", zoop_size_week12,
  "Zoop biomass",   "Composition",    zoop_comp_week12,
  "Zoop biomass",   "Phyto biomass",  zoop_phyto_week12
)


edge_summ_week12 <- edges_draws_week12 %>%
  mutate(
    Mean      = purrr::map_dbl(vec, mean),
    L89       = purrr::map_dbl(vec, ~quantile(.x, 0.055)),
    U89       = purrr::map_dbl(vec, ~quantile(.x, 0.945)),
    HDIinROPE = purrr::map_dbl(vec, ~hdi_in_rope_pct(.x, thr_week12, ci = 0.89)),
    decision  = case_when(HDIinROPE >= 95 ~ "null",
                          HDIinROPE <=  5 ~ "nonzero",TRUE            ~ "ambiguous")
    )

style_df_week12 <- edge_summ_week12 %>% # 
  transmute(from, to, Mean, L89, U89, HDIinROPE, decision,
            sign = ifelse(Mean >= 0, "pos", "neg"),
            is_sig = decision == "nonzero",
            color = case_when(is_sig & sign=="pos" ~ "#222222",
                              is_sig & sign=="neg" ~ "#C62828",
                              TRUE                 ~ "#878787"),
            style = ifelse(is_sig, "solid", "dotted"),
            penwidth = ifelse(is_sig, 5, 4),
            arrowsize = ifelse(is_sig, 1.2, 0.9))

tg12 <- dplyr::filter(style_df_week12, from == "Temp", to == "GPP") %>% dplyr::slice(1)
if (nrow(tg12) != 1) stop("Temp->GPP style not found (week 12)")

all_nodes_week12 <- unique(c(style_df_week12$from, style_df_week12$to, top_nodes, mid2_nodes, mid3_nodes, bot_nodes))
node_lines_week12 <- vapply(all_nodes_week12, function(lbl) sprintf('  "%s" [shape=ellipse];', lbl), character(1))

edge_lines_week12 <- vapply(
  seq_len(nrow(subset(style_df_week12, !(from == "Temp" & to == "GPP")))),
  function(i) {
    sdf <- subset(style_df_week12, !(from == "Temp" & to == "GPP"))
    sprintf('  "%s" -> "%s" [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s];',
            sdf$from[i], sdf$to[i], sdf$color[i], sdf$penwidth[i], sdf$arrowsize[i], sdf$style[i])
  }, character(1)
)
temp_gpp_left_week12 <- sprintf(
  '  "Temp":s -> "GPP":nw [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s, weight=70, minlen=2, constraint=false];',
  tg12$color, tg12$penwidth, tg12$arrowsize, tg12$style
)

newdot_week12 <- paste(
  dot_header,
  paste(node_lines_week12, collapse = "\n"),
  rank_top, rank_mid2_redo, rank_mid3, rank_bot,
  horiz_top, horiz_mid2, horiz_mid3, keep_N_right, center_helps,
  paste(edge_lines_week12, collapse = "\n"),
  temp_gpp_left_week12,
  paste(center_lines, collapse = "\n"),
  "}", sep = "\n"
)

SEM_fig_week12 <- DiagrammeR::grViz(newdot_week12, width = "100%")
SEM_fig_week12
svg12 <- export_svg(SEM_fig_week12)
#rsvg_pdf(charToRaw(svg12), file = "sem_diagram_week12.pdf", width = 1210, height = 590)

effect_draws_week12 <- tibble::tibble(
  Effect = c("Temp→GPP total","N→GPP total",
             "Low dispersal → GPP (ind.)","High dispersal → GPP (ind.)",
             "Size→GPP","Comp→GPP","Biomass→GPP",
             "Temp→GPP (direct)","Temp→GPP (ind.)"),
  Draws  = list(
    temp_total_week12,
    N_total_week12,
    ind_low_week12, ind_high_week12,
    size_gpp_week12, comp_gpp_week12, phyto_gpp_week12,
    temp_gpp_week12, (ind_temp_size_week12 + ind_temp_comp_week12 + ind_temp_phyto_week12)
  )
) %>% tidyr::unnest_longer(Draws, indices_include = FALSE, values_to = "value")

rope_band_week12 <- data.frame(xmin = -thr_week12, xmax = thr_week12)
lab_vec_week12 <- setNames(paste0(unique(effect_draws_week12$Effect)), unique(effect_draws_week12$Effect))

ROPEs_week12 <- ggplot(effect_draws_week12, aes(x = value)) +
  geom_rect(data = rope_band_week12,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.3, fill = "pink3") +
  geom_density(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~ Effect, scales = "free", ncol = 3,
             labeller = labeller(Effect = lab_vec_week12)) +
  labs(x = "Effect (log-GPP scale)", y = "Density",
       title = "Posterior densities with ROPE (±0.20·SD(GPP)) — Week 12") +
  theme_classic(base_size = 12)


# print(edge_summ_week4,n=22) #peeking 
# print(edge_summ_week8,n=22)
# print(edge_summ_week12,n=22)
edge_summ_all <- bind_rows(edge_summ_week4,edge_summ_week8_NMDS,edge_summ_week12_NMDS)
#print(edge_summ_all, n=63)

ROPEs_week4
ROPEs_week8
ROPEs_week12


# One ROPE figure summarizing Week-specific DAGs (in supplement)

rope_band_week4  <- data.frame(xmin = -thr_week4,  xmax =  thr_week4)
rope_band_week8  <- data.frame(xmin = -thr_week8,  xmax =  thr_week8)
rope_band_week12 <- data.frame(xmin = -thr_week12, xmax =  thr_week12)

ed4  <- effect_draws_week4  %>% mutate(Week = "Week 4")
ed8  <- effect_draws_week8  %>% mutate(Week = "Week 8")
ed12 <- effect_draws_week12 %>% mutate(Week = "Week 12")
effect_draws_all <- bind_rows(ed4, ed8, ed12)

rb4  <- rope_band_week4  %>% mutate(Week = "Week 4")
rb8  <- rope_band_week8  %>% mutate(Week = "Week 8")
rb12 <- rope_band_week12 %>% mutate(Week = "Week 12")
rope_band_all <- bind_rows(rb4, rb8, rb12)

effect_draws_all <- effect_draws_all %>%
  mutate(Week = factor(Week, levels = c("Week 4","Week 8","Week 12")))
rope_band_all <- rope_band_all %>%
  mutate(Week = factor(Week, levels = c("Week 4","Week 8","Week 12")))

ROPEs_all_20percent <- ggplot(effect_draws_all, aes(x = value)) +
  geom_rect(data = rope_band_all,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.3, fill = "pink3") +
  geom_density(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_grid(rows = vars(Week), cols = vars(Effect),
             scales = "free",  
             labeller = labeller(Effect = lab_vec_week12)) +  
  labs(x = "Effect (log-GPP scale)", y = "Density",
       title = "Posterior densities with ROPE (±0.20·SD)") +
  theme_classic(base_size = 12)+theme(strip.background = element_rect(fill = "white", color = "black"),
                                      panel.background = element_blank(),
                                      panel.border     = element_rect(color = "black", fill = NA),
                                      plot.background  = element_blank(),
                                      panel.grid       = element_blank(),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),axis.text.y=element_blank(),strip.text.y=element_text(size=12),strip.text.x=element_text(size=12))

ROPEs_all_20percent
#ggsave(ROPEs_all_20percent,file="ROPES_by_week_20percent.pdf",height=6.6,width=21)


# All weeks summary DAG (Fig. 2d) -------------------------------------------------------------
PW_CONSISTENT <- 8.0   # width when nonzero in all weeks (thick)
PW_SOMETIME   <- 2.5  # width when nonzero in 1–2 weeks (thin)
PW_NULL       <- 3.0   # width when never nonzero (thinnest)
COL_SIG       <- "#222222"  # black for any significant edge
COL_NULL      <- "#878787"  # grey for null
ARROWSIZE_ALL <- 1.2        # uniform arrowhead size for all edges

# Combine week summaries 
edges_allweeks <- bind_rows(
  mutate(edge_summ_week4,  week = "4"),
  mutate(edge_summ_week8,  week = "8"),
  mutate(edge_summ_week12, week = "12")
)
print(edges_allweeks,n=63 )

synth_map <- edges_allweeks %>%
  group_by(from, to) %>%
  summarise(n_nonzero = sum(decision == "nonzero", na.rm = TRUE), .groups = "drop") %>%
  mutate(
    cat = case_when(
      n_nonzero == 3 ~ "consistent",
      n_nonzero >= 1 ~ "sometimes",
      TRUE           ~ "insig"
    ),
    color = ifelse(cat == "insig", COL_NULL, COL_SIG),
    style = ifelse(cat == "insig", "dotted", "solid"),
    penwidth = case_when(
      cat == "consistent" ~ PW_CONSISTENT,
      cat == "sometimes"  ~ PW_SOMETIME,
      TRUE                ~ PW_NULL
    ),
    arrowsize = ARROWSIZE_ALL
  ) %>%
  select(from, to, color, style, penwidth, arrowsize)

# Take working week-4 style_df and override color/style/penwidth/arrowsize
style_df_all <- style_df_week4 %>%              
  select(from, to) %>% distinct() %>%           
  left_join(synth_map, by = c("from","to"))

tg4 <- dplyr::filter(style_df_all, from == "Temp", to == "GPP") %>% dplyr::slice(1)
if (nrow(tg4) != 1) stop("Temp->GPP style not found (synth)")

top_nodes <- c("Temp", "Dispersal: low", "Dispersal: high")
mid2_nodes <- c("Zoop biomass", "Total N")
mid3_nodes <- c("Phyto biomass", "Composition", "Size structure")
bot_nodes <- c("GPP")

dot_header <- paste(
  'digraph sem {',
  '  graph [layout=dot, rankdir=TB, newrank=true,',
  '         splines=curved, overlap=false,',
  '         nodesep="0.8", ranksep="1.8",',
  '         margin=0.05, pad="0.1,0.1", center=true];',
  '  node  [shape=ellipse, style="filled", fillcolor="white", color="black",',
  '         fontname="Helvetica", fontsize=20, fixedsize=true, width=2, height=1];',
  '  edge  [fontname="Helvetica", fontsize=14, color="#888888", arrowsize=0.9, labelfloat=false];',
  sep = "\n"
)

all_nodes_allweeks <- unique(c(style_df_all$from, style_df_all$to, top_nodes, mid2_nodes, mid3_nodes, bot_nodes))
node_lines_allweeks <- vapply(all_nodes_allweeks, function(lbl) sprintf('  "%s" [shape=ellipse];', lbl), character(1))

rank_top <- paste(
  '  { rank=same;',
  '    Ctop [shape=point, width=0, label="", style=invis];',
  '    SpacerTopL [shape=point, label="", style=invis, width=0.6, height=0.01];',
  '    "Temp";',
  '    "Dispersal: low";',
  '    "Dispersal: high";',
  '  }', sep = "\n"
)
rank_mid2_redo <- paste(
  '  { rank=same;',
  '    Cmid2 [shape=point, width=0, label="", style=invis];',
  '    "Zoop biomass";',
  '    SpacerN [shape=point, label="", style=invis, width=0.6, height=0.01];',
  '    "Total N";',
  '  }', sep = "\n"
)
rank_mid3 <- paste(
  '  { rank=same;',
  '    Cmid3 [shape=point, width=0, label="", style=invis];',
  '    "Size structure";',
  '    Spacer3a [shape=point, label="", style=invis, width=0.6, height=0.01];',
  '    "Composition";',
  '    Spacer3b [shape=point, label="", style=invis, width=0.6, height=0.01];',
  '    "Phyto biomass";',
  '  }', sep = "\n"
)
rank_bot <- paste(
  '  { rank=same;',
  '    Cbot [shape=point, width=0, label="", style=invis];',
  '    "GPP";',
  '  }', sep = "\n"
)
horiz_top <- paste(
  '  { rank=same;',
  '    SpacerTopL -> "Temp" -> "Dispersal: low" -> "Dispersal: high" [style=invis, weight=80];',
  '  }', sep = "\n"
)
horiz_mid2 <- paste(
  '  { rank=same;',
  '    "Zoop biomass" -> SpacerN -> "Total N" [style=invis, weight=70];',
  '  }', sep = "\n"
)
horiz_mid3 <- paste(
  '  { rank=same;',
  '    "Size structure" -> Spacer3a -> "Composition" -> Spacer3b -> "Phyto biomass" [style=invis, weight=85];',
  '  }', sep = "\n"
)
keep_N_right <- paste(
  '  "Dispersal: high" -> "Total N" [style=invis, weight=51];',
  '  "Phyto biomass"   -> "Total N" [style=invis, weight=47];', sep = "\n"
)
center_helps <- paste(
  '  "Dispersal: low"  -> "Zoop biomass" [style=invis, weight=20];',
  '  "Dispersal: high" -> "Zoop biomass" [style=invis, weight=20];',
  '  "Size structure"  -> "Zoop biomass" [style=invis, weight=20];',
  '  "Phyto biomass"   -> "Zoop biomass" [style=invis, weight=20];', sep = "\n"
)
center_lines <- c(
  '  Ctop -> Cmid2 -> Cmid3 -> Cbot [style=invis, weight=200];',
  paste(sprintf('  Ctop  -> "%s" [style=invis, constraint=false, weight=80];', top_nodes), collapse = "\n"),
  paste(sprintf('  Cmid2 -> "%s" [style=invis, constraint=false, weight=80];', mid2_nodes), collapse = "\n"),
  paste(sprintf('  Cmid3 -> "%s" [style=invis, constraint=false, weight=80];', mid3_nodes), collapse = "\n"),
  '  Cbot  -> "GPP" [style=invis, constraint=false, weight=80];'
)

style_df_no_TG_allweeks <- subset(style_df_all, !(from == "Temp" & to == "GPP"))
edge_lines_allweeks <- vapply(
  seq_len(nrow(style_df_no_TG_allweeks)),
  function(i) sprintf(
    '  "%s" -> "%s" [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s];',
    style_df_no_TG_week4$from[i], style_df_no_TG_week4$to[i],
    style_df_no_TG_week4$color[i], style_df_no_TG_week4$penwidth[i],
    style_df_no_TG_week4$arrowsize[i], style_df_no_TG_week4$style[i]
  ), character(1)
)
temp_gpp_left_allweeks <- sprintf(
  '  "Temp":s -> "GPP":nw [color="%s", penwidth=%.2f, arrowsize=%.2f, style=%s, weight=70, minlen=2, constraint=false];',
  tg4$color, tg4$penwidth, tg4$arrowsize, tg4$style
)

newdot_allweeks <- paste(
  dot_header,
  paste(node_lines_allweeks, collapse = "\n"),
  rank_top, rank_mid2_redo, rank_mid3, rank_bot,
  horiz_top, horiz_mid2, horiz_mid3, keep_N_right, center_helps,
  paste(edge_lines_allweeks, collapse = "\n"),
  temp_gpp_left_allweeks,
  paste(center_lines, collapse = "\n"),
  "}", sep = "\n"
)

SEM_fig_ALLWEEKS_summary <- DiagrammeR::grViz(newdot_allweeks, width = "100%")
SEM_fig_ALLWEEKS_summary

svgSEM_fig_ALLWEEKS_summary_20percent <- export_svg(SEM_fig_ALLWEEKS_summary)
#rsvg_pdf(charToRaw(svgSEM_fig_ALLWEEKS_summary_20percent), file = "SEM_DAG_allweeks_20percent.pdf", width = 700, height = 590)


-----------------------------------------------------
# Checking sensitivity to ROPE width (5%,10%,20%) for Table S3

rope_sensitivity_df <- function(edge_df, sdy) {
  widths <- c(0.05, 0.10, 0.20)
  purrr::map_dfr(
    seq_len(nrow(edge_df)), \(i) {
      draws <- edge_df$vec[[i]]
      purrr::map_dfc(widths, \(w) {
        thr <- w * sdy
        pct <- hdi_in_rope_pct(draws, thr, ci = 0.89)
        tibble(!!paste0("ROPE_", w*100, "_pct") := pct)
      }) |>
        mutate(Effect = paste(edge_df$from[i], "→", edge_df$to[i]))
    }
  )
}

sdy_week4  <- sd(df_sem_subset_week4$ln_GPP_umol_L_h,  na.rm = TRUE)
sdy_week8  <- sd(df_sem_subset_week8$ln_GPP_umol_L_h,  na.rm = TRUE)
sdy_week12 <- sd(df_sem_subset_week12$ln_GPP_umol_L_h, na.rm = TRUE)

rope_sens_week4  <- rope_sensitivity_df(edge_summ_week4,  sdy_week4)
rope_sens_week8  <- rope_sensitivity_df(edge_summ_week8_NMDS,  sdy_week8)
rope_sens_week12 <- rope_sensitivity_df(edge_summ_week12_NMDS, sdy_week12)

rope_sens_week4  <- rope_sens_week4  |> mutate(Week = "4")
rope_sens_week8  <- rope_sens_week8  |> mutate(Week = "8")
rope_sens_week12 <- rope_sens_week12 |> mutate(Week = "12")

rope_sensitivity_allweeks <- bind_rows(rope_sens_week4, rope_sens_week8, rope_sens_week12)

rope_sensitivity_allweeks <- rope_sensitivity_allweeks |>
  mutate(
    sensitive_to_rope = (ROPE_5_pct < 5 & (ROPE_10_pct >= 5 | ROPE_20_pct >= 5)) |
      (ROPE_10_pct < 5 & ROPE_20_pct >= 5)
  )
