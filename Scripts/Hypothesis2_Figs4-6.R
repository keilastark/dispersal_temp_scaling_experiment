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

exp_data<-read.csv("metacom_exp_data_final.csv")
# make sure factors properly set
exp_data$metacommunity <- factor(exp_data$metacommunity)
exp_data$wattage <- factor(exp_data$wattage, levels=c("0","100","200","300"))
exp_data$dispersal <- factor(exp_data$dispersal, levels=c("none","low","high"))
exp_data$tank <- factor(exp_data$tank)
exp_data$sample_date <- as.character(exp_data$sample_date)

############################################################
# HYPOTHESIS 2A: BIOMASS ------------------------------------
############################################################

#Manuscript Fig 4:
# Time series of GPP and biomass by wattage and dispersal treatments:

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
  ) + scale_y_continuous(limits=c(1.5,3) 
  )+
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

biom_v_time <- exp_data %>%
  ggplot(aes(x = week_number, y = ln_fluoro_ug_C_L, color = wattage, group = wattage, shape = wattage)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 0.6) +
  facet_wrap(~dispersal, labeller = labeller(dispersal = dispersal_labels)) +
  scale_color_brewer(name = "Wattage treatment", palette = "RdBu", direction = -1) +
  scale_shape_manual(name = "Wattage treatment",values = c(16, 17, 15, 18)) + 
  theme_classic() +
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

legend_gpp <- ggpubr::get_legend(
  gpp_v_time + theme(legend.position = "bottom")
)


fig4main <- plot_grid(
  gpp_nolegend,biomass_nolegend,
  labels=c("(a)","(b)"),label_size = 15, vjust=-0.5,hjust=1.3,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1)
)

fig4<- plot_grid(
  fig4main,legend_gpp,labels=NULL,
  ncol = 1,
  rel_heights = c(1,0.1)
)+
  theme(plot.margin = margin(t = 25, r = 0, b = 0, l = 30))

#ggsave(fig4,file="biom_gpp_timeseries.pdf",width=7.5,height=6.6)


################# Levene's tests of residual variance by dispersal group
# Seeing if dispersal is associated with differences in variance in GPP and biomass 
# uses car package

w2 <- exp_data %>% filter(week_number == "2")
w4 <- exp_data %>% filter(week_number == "4")
w6 <- exp_data %>% filter(week_number == "6")
w8 <- exp_data %>% filter(week_number == "8")
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
fit8_gpp <- lm(ln_GPP_umol_L_h  ~ dispersal, data = w8)
leveneTest(fit8_bio, center = mean)
leveneTest(fit8_gpp, center = mean)

fit10_bio <- lm(ln_fluoro_ug_C_L ~ dispersal, data = w10)
fit10_gpp <- lm(ln_GPP_umol_L_h  ~ dispersal, data = w10)
leveneTest(fit10_bio,   center = mean)
leveneTest(fit10_gpp,  center = mean)

fit12_bio <- lm(ln_fluoro_ug_C_L ~ dispersal, data = w12)
fit12_gpp <- lm(ln_GPP_umol_L_h  ~ dispersal, data = w12)
leveneTest(fit12_bio, center = mean)
leveneTest(fit12_gpp, center = mean)

pairwise_levene <- function(formula, data, center = mean, adjust = "holm") {
  resp <- all.vars(formula)[1]
  grp  <- all.vars(formula)[2]
  levs <- levels(droplevels(as.factor(data[[grp]])))
  pr   <- combn(levs, 2, simplify = FALSE)
  
  out <- map_dfr(pr, function(p) {
    sub <- droplevels(data[data[[grp]] %in% p, , drop = FALSE])
    a   <- car::leveneTest(formula, data = sub, center = center)
    
    sd1 <- sd(sub[[resp]][sub[[grp]] == p[1]], na.rm = TRUE)
    sd2 <- sd(sub[[resp]][sub[[grp]] == p[2]], na.rm = TRUE)
    
    tibble(
      g1 = p[1], g2 = p[2],
      df1 = a$Df[1], df2 = a$Df[2],
      F   = a$`F value`[1],
      p   = a$`Pr(>F)`[1],
      n1  = sum(sub[[grp]] == p[1]),
      n2  = sum(sub[[grp]] == p[2]),
      sd1 = sd1, sd2 = sd2,
      var_ratio = (sd1^2) / (sd2^2)
    )
  }) %>% mutate(p_adj = p.adjust(p, method = adjust))
  
  out
}

weeks <- list(w2 = w2, w4 = w4, w6 = w6, w8 = w8, w10 = w10, w12 = w12)

bio_pw <- imap(weeks, ~ pairwise_levene(ln_fluoro_ug_C_L ~ dispersal, .x, center = mean)) %>%
  bind_rows(.id = "week")

gpp_pw <- imap(weeks, ~ pairwise_levene(ln_GPP_umol_L_h ~ dispersal, .x, center = mean)) %>%
  bind_rows(.id = "week")

bio_pw
gpp_pw
# can use these outputs to color non-sig points in fig.3 grey

############################################################
# HYPOTHESIS 2B: PHYTOPLANKTON COMMUNITY SIZE STRUCTURE -------------------
############################################################

## Generating frequency plot of size bins (Fig. 5a)

size_cols <- c("X0_5um","X5_10um","X10_15um","X15_20um",
               "X20_25um","X25_30um","X30_35um","X35_40um")

phyto_size_na_rm <- exp_data %>%
  filter(!if_all(all_of(size_cols), is.na))

phyto_size_struc <- phyto_size_na_rm %>%
  dplyr::select(c(all_of(size_cols),wattage,dispersal,week_number)) 

phyto_size_long <- phyto_size_struc %>%
  mutate(row_sum = rowSums(across(starts_with("X")))) %>%
  mutate(across(starts_with("X"), ~ .x / row_sum)) %>%
  select(-row_sum) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "size_bin",
    values_to = "prop"
  ) %>%
  mutate(
    size_bin = gsub("^X", "", size_bin),
    size_bin = gsub("um", " µm", size_bin),
    size_bin = factor(size_bin, levels = c(
      "0_5 µm", "5_10 µm", "10_15 µm", "15_20 µm",
      "20_25 µm", "25_30 µm", "30_35 µm", "35_40 µm"
    ))
  )
rdbu_colors <- brewer.pal(4, "RdBu") %>% rev()
names(rdbu_colors) <- c("0", "100", "200", "300")

phyto_summary_watt <- phyto_size_long %>%
  group_by(week_number, wattage, size_bin) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE),
            se_prop   = sd(prop, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

size_breaks <- c("0_5 µm","5_10 µm","10_15 µm","15_20 µm",
                 "20_25 µm","25_30 µm","30_35 µm","35_40 µm")

size_labels_parse <- c("0<=5","5<=10","10<=15","15<=20",
                       "20<=25","25<=30","30<=35","35<=40")
lab_vec <- setNames(size_labels_parse, size_breaks)

# Size histogram by wattage (not in MS)
freq_wattage <- ggplot(phyto_summary_watt, aes(x = size_bin, y = mean_prop, fill = wattage)) +
  geom_col(position = position_dodge(width = 0.9), color = "black", width = 0.9) +
  geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~ week_number, nrow = 2,
             labeller = labeller(week_number = c("2"="Week 2","4"="Week 4","6"="Week 6",
                                                 "8"="Week 8","10"="Week 10","12"="Week 12"))) +
  scale_x_discrete(breaks = size_breaks,
                   labels = function(x) parse(text = lab_vec[x])) +
  
  scale_fill_brewer(palette="RdBu",direction=-1, name = "Wattage") +
  labs(
    x = "Phytoplankton size bin (ESD, µm)",
    y = "Mean proportion (±SE)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 13, family = "Helvetica"),
    strip.text        = element_text(size = 12, family = "Helvetica"),
    legend.text = element_text(size = 12, family = "Helvetica"),
    axis.title = element_text(size = 14, family = "Helvetica"),
    axis.text = element_text(size = 12, family = "Helvetica"),
    axis.title.x = element_text(vjust=-1),
    axis.title.y = element_text(vjust=2),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  )
# Supplementary Fig: 
#ggsave(freq_wattage,file="watt_freq_size_plot.pdf",width=6.1,height=4.8)


#### Now dispersal size histogram for manuscript ###

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

#size_breaks <- c("0_5 µm","5_10 µm","10_15 µm","15_20 µm", "20_25 µm","25_30 µm","30_35 µm","35_40 µm")
#size_labels_parse <- c("0<=5","5<=10","10<=15","15<=20", "20<=25","25<=30","30<=35","35<=40")
#lab_vec <- setNames(size_labels_parse, size_breaks)

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


### Manuscript figure 6
#ggsave(disp_freq_size_plot,file="disp_freq_size_plot.pdf",width=6.1,height=4.8)

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

combined_size <- disp_freq_size_plot  + plot_spacer()+ phyto_size_watt_time + 
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

ggsave(combined_size, file="size_figs.pdf",width=12.9,height=4.6)

# Supplementary zoop biomass vs inedible phyto abundance figure

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


###### Multivariate analyses for size structure and taxonomic composition --------

size_cols <- c("X0_5um","X5_10um","X10_15um","X15_20um",
               "X20_25um","X25_30um","X30_35um","X35_40um")

phyto_size_na_rm <- exp_data %>%
  filter(!if_all(all_of(size_cols), is.na))

phyto_size_struc_env <- phyto_size_na_rm %>%
  select(sample_date, tank, wattage, temp, week_number,dispersal,metacommunity) 

phyto_size_env$wattage<-factor(phyto_size_env$wattage, levels=c("0","100","200","300"))


phyto_size_struc_mat <- phyto_size_na_rm %>%
  dplyr::select(c(all_of(size_cols))) 

row_key <- paste(phyto_size_na_rm$week_number, phyto_size_na_rm$tank, sep = "__")
rownames(phyto_size_struc_mat) <- row_key

phyto_size_nmds_result <- metaMDS(phyto_size_struc_mat, distance = "bray")
phyto_size_scores <- scores(phyto_size_nmds_result, display = "sites")
phyto_size_scores <- as.data.frame(phyto_size_scores)

# combine with the full dataset (ie for SEM analyses)

#size_scores <- phyto_size_nmds_scores %>%
#  mutate(sample_date = as.character(sample_date),
#         tank        = as.character(tank)) %>%
#  distinct(sample_date, tank, .keep_all = TRUE) %>%
#  transmute(sample_date, tank,
#            size_NMDS1 = NMDS1,
#            size_NMDS2 = NMDS2)
#
#comp_scores <- phyto_comp_nmds_scores %>%
#  mutate(sample_date = as.character(sample_date),
#         tank        = as.character(tank)) %>%
#  distinct(sample_date, tank, .keep_all = TRUE) %>%
#  transmute(sample_date, tank,
#            comp_NMDS1 = NMDS1,
#            comp_NMDS2 = NMDS2)
#
#exp_data <- exp_data %>%
#  left_join(size_scores, by = c("sample_date", "tank")) %>%
#  left_join(comp_scores, by = c("sample_date", "tank"))


### Size PERMANOVA, all weeks together

phyto_size_perman <- adonis2(
  phyto_size_struc_mat ~ week_number + temp +dispersal,
  by = "terms",
  strata = phyto_size_struc_env$metacommunity,
  method = "bray",
  data = phyto_size_struc_env
)

### Size PERMANOVAs by week

### Filter week-specific data
week2phyto_size_struc_mat <-   phyto_size_na_rm %>%
  filter(week_number==2) %>%
  select(X0_5um, X5_10um, X10_15um, X15_20um, 
         X20_25um, X25_30um, X30_35um, X35_40um) 

week2phyto_size_struc_env <-   phyto_size_na_rm %>%
  filter(week_number==2) %>%
  select(tank,wattage,metacommunity,dispersal,temp,sample_date,week_number) 

week4phyto_size_struc_mat <- phyto_size_na_rm %>%
  filter(week_number == 4) %>%
  select(X0_5um, X5_10um, X10_15um, X15_20um, 
         X20_25um, X25_30um, X30_35um, X35_40um) 

week4phyto_size_struc_env <- phyto_size_na_rm %>%
  filter(week_number == 4) %>%
  select(tank, wattage, metacommunity, dispersal, temp, sample_date, week_number)

week6phyto_size_struc_mat <- phyto_size_na_rm %>%
  filter(week_number == 6) %>%
  select(X0_5um, X5_10um, X10_15um, X15_20um, 
         X20_25um, X25_30um, X30_35um, X35_40um) 

week6phyto_size_struc_env <- phyto_size_na_rm %>%
  filter(week_number == 6) %>%
  select(tank, wattage, metacommunity, dispersal, temp, sample_date, week_number)

week8phyto_size_struc_mat <- phyto_size_na_rm %>%
  filter(week_number == 8) %>%
  select(X0_5um, X5_10um, X10_15um, X15_20um, 
         X20_25um, X25_30um, X30_35um, X35_40um) 

week8phyto_size_struc_env <- phyto_size_na_rm %>%
  filter(week_number == 8) %>%
  select(tank, wattage, metacommunity, dispersal, temp, sample_date, week_number)

week10phyto_size_struc_mat <- phyto_size_na_rm %>%
  filter(week_number == 10) %>%
  select(X0_5um, X5_10um, X10_15um, X15_20um, 
         X20_25um, X25_30um, X30_35um, X35_40um) 

week10phyto_size_struc_env <- phyto_size_na_rm %>%
  filter(week_number == 10) %>%
  select(tank, wattage, metacommunity, dispersal, temp, sample_date, week_number)

week12phyto_size_struc_mat <- phyto_size_na_rm %>%
  filter(week_number == 12) %>%
  select(X0_5um, X5_10um, X10_15um, X15_20um, 
         X20_25um, X25_30um, X30_35um, X35_40um) 

week12phyto_size_struc_env <- phyto_size_na_rm %>%
  filter(week_number == 12) %>%
  select(tank, wattage, metacommunity, dispersal, temp, sample_date, week_number)

### run multivariate analyses
week2_size_permanova <- adonis2(
  week2phyto_size_struc_mat ~ temp * dispersal,
  by = "terms",
  strata = week2phyto_size_struc_env$metacommunity,
  method = "bray",
  data = week2phyto_size_struc_env
)
#Df SumOfSqs      R2      F Pr(>F)  
#temp            1   0.3250 0.08017 3.4099  0.031 *
#  dispersal       2   0.2094 0.05165 1.0985  0.967  
#temp:dispersal  2   0.0882 0.02176 0.4628  0.825  
#Residual       36   3.4314 0.84641                
#Total          41   4.0540 1.00000    

week4_size_permanova <- adonis2(
  week4phyto_size_struc_mat ~ temp * dispersal,
  by = "terms",
  strata = week4phyto_size_struc_env$metacommunity,
  method = "bray",
  data = week4phyto_size_struc_env
)
#Df SumOfSqs      R2      F Pr(>F)   
#temp            1   0.3340 0.08540 4.4535  0.008 **
#  dispersal       2   0.3896 0.09960 2.5969  0.265   
#temp:dispersal  2   0.1125 0.02877 0.7500  0.567   
#Residual       41   3.0753 0.78623                 
#Total          46   3.9114 1.00000    

week6_size_permanova <- adonis2(
  week6phyto_size_struc_mat ~ temp * dispersal,
  by = "terms",
  strata = week6phyto_size_struc_env$metacommunity,
  method = "bray",
  data = week6phyto_size_struc_env
)
#Df SumOfSqs      R2      F Pr(>F)  
#temp            1   0.1804 0.05557 2.8304  0.023 *
#  dispersal       2   0.2759 0.08499 2.1642  0.014 *
#  temp:dispersal  2   0.1129 0.03478 0.8858  0.370  
#Residual       42   2.6773 0.82466                
#Total          47   3.2466 1.00000   

week8_size_permanova <- adonis2(
  week8phyto_size_struc_mat ~ temp * dispersal,
  by = "terms",
  strata = week8phyto_size_struc_env$metacommunity,
  method = "bray",
  data = week8phyto_size_struc_env
)
#Df SumOfSqs      R2      F Pr(>F)
#temp            1   0.0833 0.02349 1.2303  0.262
#dispersal       2   0.7182 0.20255 5.3054  0.727
#temp:dispersal  2   0.1046 0.02949 0.7723  0.500
#Residual       39   2.6398 0.74448              
#Total          44   3.5459 1.00000          

# Week 10
week10_size_permanova <- adonis2(
  week10phyto_size_struc_mat ~ temp * dispersal,
  by = "terms",
  strata = week10phyto_size_struc_env$metacommunity,
  method = "bray",
  data = week10phyto_size_struc_env
)
#Df SumOfSqs      R2      F Pr(>F)   
#temp            1   0.1164 0.03316 1.8393  0.137   
#dispersal       2   0.3877 0.11044 3.0628  0.010 **
#  temp:dispersal  2   0.4116 0.11724 3.2514  0.012 * 
#  Residual       41   2.5951 0.73917                 
#Total          46   3.5109 1.00000      

week12_size_permanova <- adonis2(
  week12phyto_size_struc_mat ~ temp * dispersal,
  by = "terms",
  strata = week12phyto_size_struc_env$metacommunity,
  method = "bray",
  data = week12phyto_size_struc_env
)

#Df SumOfSqs      R2      F Pr(>F)
#temp            1   0.1241 0.03281 1.5204  0.170
#dispersal       2   0.4760 0.12590 2.9171  0.816
#temp:dispersal  2   0.0805 0.02128 0.4931  0.739
#Residual       38   3.1006 0.82001              
#Total          43   3.7811 1.00000       


###### PERMDISPs: Community size structure --------------------------------------

## PERMDISP by week number

# week 2
phyto_size_dist_mat2 <- vegdist(week2phyto_size_struc_mat, method = "bray")

phyto_size_bd2 <- betadisper(phyto_size_dist_mat2, type="centroid",group = week2phyto_size_struc_env$metacommunity) 

phyto_size_bd_df2 <- data.frame(
  tank      = week2phyto_size_struc_env$tank,
  metacom         = phyto_size_bd2$group, 
  dispersal= week2phyto_size_struc_env$dispersal,
  DistToCentroid = phyto_size_bd2$distances # distance to centroid of metacom
)

phyto_week2_dispersion <- aov(DistToCentroid ~ dispersal, data = phyto_size_bd_df2)
anova(phyto_week2_dispersion)
#Analysis of Variance Table
#
#Response: DistToCentroid
#Df  Sum Sq  Mean Sq F value Pr(>F)
#dispersal  2 0.02765 0.013825  0.8266 0.4445
#Residuals 42 0.70247 0.016726   

# week 4
phyto_size_dist_mat4 <- vegdist(week4phyto_size_struc_mat, method = "bray")

phyto_size_bd4 <- betadisper(phyto_size_dist_mat4, type="centroid",group = week4phyto_size_struc_env$metacommunity) 

phyto_size_bd_df4 <- data.frame(
  tank      = week4phyto_size_struc_env$tank,
  metacom         = phyto_size_bd4$group, 
  dispersal= week4phyto_size_struc_env$dispersal,
  DistToCentroid = phyto_size_bd4$distances
)

phyto_week4_dispersion <- aov(DistToCentroid ~ dispersal, data = phyto_size_bd_df4)
anova(phyto_week4_dispersion) 

#Analysis of Variance Table
#
#Response: DistToCentroid
#Df  Sum Sq  Mean Sq F value   Pr(>F)   
#dispersal  2 0.17119 0.085597  7.5362 0.001533 **
#  Residuals 44 0.49976 0.011358                    


# week 6 
phyto_size_dist_mat6<- vegdist(week6phyto_size_struc_mat, method = "bray")

phyto_size_bd6 <- betadisper(phyto_size_dist_mat6, group = week6phyto_size_struc_env$metacommunity,type="centroid")

phyto_size_bd_df6 <- data.frame(
  tank      = week6phyto_size_struc_env$tank,
  metacom         = phyto_size_bd6$group, 
  dispersal= week6phyto_size_struc_env$dispersal,
  DistToCentroid = phyto_size_bd6$distances
)

phyto_week6_dispersion <- aov(DistToCentroid ~ dispersal, data = phyto_size_bd_df6)
anova(phyto_week6_dispersion)

#Analysis of Variance Table
#
#Response: DistToCentroid
#Df  Sum Sq   Mean Sq F value Pr(>F)
#dispersal  2 0.04540 0.0227024  2.3846 0.1037
#Residuals 45 0.42841 0.0095203         


# week 8

phyto_size_dist_mat8 <- vegdist(week8phyto_size_struc_mat, method = "bray")

phyto_size_bd8 <- betadisper(phyto_size_dist_mat8, group = week8phyto_size_struc_env$metacommunity,type="centroid")

phyto_size_bd_df8 <- data.frame(
  tank      = week8phyto_size_struc_env$tank,
  metacom         = phyto_size_bd8$group, 
  dispersal= week8phyto_size_struc_env$dispersal,
  DistToCentroid = phyto_size_bd8$distances
)

phyto_week8_dispersion <- aov(DistToCentroid ~ dispersal, data = phyto_size_bd_df8)
anova(phyto_week8_dispersion) 

#Response: DistToCentroid
#Df   Sum Sq   Mean Sq F value  Pr(>F)  
#dispersal  2 0.058714 0.0293571   5.078 0.01059 *


# week 10 

phyto_size_dist_mat10 <- vegdist(week10phyto_size_struc_mat, method = "bray")

phyto_size_bd10 <- betadisper(phyto_size_dist_mat10, group = week10phyto_size_struc_env$metacommunity,type="centroid")

phyto_size_bd_df10 <- data.frame(
  tank      = week10phyto_size_struc_env$tank,
  temp  = week10phyto_size_struc_env$temp,
  metacom         = phyto_size_bd10$group, 
  dispersal= week10phyto_size_struc_env$dispersal,
  DistToCentroid = phyto_size_bd10$distances
)

phyto_week10_dispersion <- aov(DistToCentroid ~ dispersal, data = phyto_size_bd_df10)
anova(phyto_week10_dispersion) 
#Analysis of Variance Table
#
#Response: DistToCentroid
#Df  Sum Sq  Mean Sq F value  Pr(>F)  
#dispersal  2 0.10115 0.050574  4.0579 0.02413 *


# week 12 

phyto_size_dist_mat12 <- vegdist(week12phyto_size_struc_mat, method = "bray")

phyto_size_bd12 <- betadisper(phyto_size_dist_mat12, group = week12phyto_size_struc_env$metacommunity,type="centroid")

phyto_size_bd_df12 <- data.frame(
  tank      = week12phyto_size_struc_env$tank,
  temp  = week12phyto_size_struc_env$temp,
  metacom         = phyto_size_bd12$group, 
  dispersal= week12phyto_size_struc_env$dispersal,
  DistToCentroid = phyto_size_bd12$distances
)

phyto_week12_dispersion <- aov(DistToCentroid ~ dispersal, data = phyto_size_bd_df12)
anova(phyto_week12_dispersion) 
#Analysis of Variance Table
#
#Response: DistToCentroid
#Df  Sum Sq  Mean Sq F value Pr(>F)
#dispersal  2 0.03554 0.017771  0.7249 0.4905
#Residuals 41 1.00513 0.024515     


###### Combining all phytoplankton size structure PERMDISPS for manuscript fig6

phyto_size_bd_all <- bind_rows(
  mutate(phyto_size_bd_df2 , week = 2),
  mutate(phyto_size_bd_df4 , week = 4),
  mutate(phyto_size_bd_df6 , week = 6),
  mutate(phyto_size_bd_df8 , week = 8),
  mutate(phyto_size_bd_df10 , week = 10),
  mutate(phyto_size_bd_df12, week = 12)
) %>% 
  relocate(week, .after = tank)  

phyto_size_bd_all <- phyto_size_bd_all %>%
  mutate(dispersal = factor(dispersal, levels = c("none", "low", "high")))

phyto_permdisp_size_all <-
  ggplot(phyto_size_bd_all,
         aes(x = factor(week),   
             y = DistToCentroid,
             group = interaction(week, dispersal),
             fill  = dispersal)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("#003f5c", "#ff6361", "#ffa600")) +
  scale_colour_manual(values = c("#003f5c", "#ff6361", "#ffa600")) +
  scale_x_discrete(name = "Week", labels = c("2","4","6", "8","10", "12")) +
  labs(title="Phytoplankton community size structure",
    y = "Distance to metacommunity centroid (PERMDISP)"
  ) +
  
  theme_classic()+
  theme(
    legend.position = "bottom",
    panel.border    = element_rect(colour = "black", fill = NA),
    axis.title      = element_text(size = 13),
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12),
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 13) 
  )+ylim(0.02,0.65)


#ggsave(phyto_permdisp_size_all, file="size_struc_permdisp.png",width=5.3,height=5.1)

# Significance codes (i.e., a, ab, b) from post-hoc comparisons
letters_df <- phyto_size_bd_all %>% 
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

letters_df <- letters_df %>% 
  group_by(week) %>% 
  mutate(letter = if (all(letter == letter[1])) "a" else letter) %>% 
  ungroup()

text_pos <- phyto_size_bd_all %>% 
  group_by(week, dispersal) %>% 
  summarise(y_height = max(DistToCentroid) + 0.02, .groups = "drop")

letters_df <- letters_df %>% left_join(text_pos, by = c("week", "dispersal"))
letters_df <- letters_df %>%           
  left_join(
    phyto_size_bd_all %>%                    
      group_by(week, dispersal) %>%
      summarise(y = max(DistToCentroid, na.rm = TRUE) + 0.02,
                .groups = "drop"),
    by = c("week", "dispersal")
  )

phyto_permdisp_size_all +
  geom_text(data = letters_df,
            aes(x = factor(week), y = y_height, label = letter, group = dispersal),
            position = position_dodge(width = 0.7),
            size = 5,
            vjust = 0)  

############################################################
# HYPOTHESIS 2C: TAXONOMIC COMPOSITION ------------------------------------
############################################################

# Read in and wrangle composition data for all three trophic levels------------------------------------------------

# Read/clean phytoplankton count (flowcam) data

Phyto_data<-read.csv("phytoplankton_diversity_data.csv")

metadata <- exp_data %>%
  mutate(
    sample_date = as.Date(sample_date),                      
    across(c(wattage, metacommunity, dispersal, temp), ~ .)
  ) %>%
  select(sample_date, tank, wattage, metacommunity, dispersal, temp)

metadata$sample_date <- factor(metadata$sample_date)
Phyto_data$tank <- factor(Phyto_data$tank)

Phyto_combined <- Phyto_data %>%
  left_join(metadata, by = c("sample_date", "tank"))

sp_cols <- which(names(Phyto_combined) == "Ceratium") :
  which(names(Phyto_combined) == "Cryptomonas.sp.")

Phyto_combined <- Phyto_combined %>%
  tidyr::drop_na( all_of(sp_cols) ) 

Phyto_combined <- Phyto_combined %>%
  mutate(week_number = dplyr::recode(sample_date, !!!week_mapping)) 
Phyto_combined$week_number <- factor(Phyto_combined$week_number)

Phyto_combined$dispersal <- factor(Phyto_combined$dispersal, ordered = TRUE)
Phyto_combined$tank <- factor(Phyto_combined$tank, ordered = TRUE)

phyto_mat <- data.matrix(Phyto_combined %>% dplyr::select(Ceratium:Cryptomonas.sp.))

env_phyto <- Phyto_combined %>%
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


#### Read zoop data
Zoop_data <- read.csv("Zoop_data.csv")

#### Read /clean bacteria data
Bacteria_dat<-read.csv("Bacterial_16S_data.csv")
Bacteria_dat$dispersal <- factor(Bacteria_dat$dispersal, ordered = TRUE)
Bacteria_dat$tank <- as.factor(Bacteria_dat$tank)


# ZOOPLANKTON -------------------------------------------------------------

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

# Run zooplankton PERMANOVAs (same result as in Thompson et al 2024 preprint)

zoop_perman4 <-adonis2(week4zoop ~dispersal*temp_mean,  method = "bray", strata = week4zoop_env$metacommunity,  by="terms", data = week4zoop_env)

zoop_perman8 <-adonis2(week8zoop ~dispersal*temp_mean,  method = "bray", strata = week8zoop_env$metacommunity,  by="terms", data = week8zoop_env)

zoop_perman12 <-adonis2(week12zoop ~dispersal*temp_mean,  method = "bray", strata = week12zoop_env$metacommunity,  by="terms", data = week12zoop_env)


### PERMDISP

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

bd_df_combined <- rbind(bd_df4,bd_df8, bd_df12)
bd_df_combined$week <- factor(bd_df_combined$week, levels = c("Week 4", "Week 8", "Week 12"))
bd_df_combined$dispersal <- factor(bd_df_combined$dispersal, levels = c("none", "low", "high"))

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


# PHYTOPLANKTON COMPOSITION -----------------------------------------------

# separate by week
week4phyto <- Phyto_combined%>%
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
phyto_perman4 <-adonis2(decostand(week4phyto,method="hellinger") ~dispersal*temp,  method = "bray", strata = week4phyto_env$metacommunity,  by="terms", data = week4phyto_env)

phyto_perman8 <-adonis2(decostand(week8phyto,method="hellinger") ~dispersal*temp,  method = "bray", strata = week8phyto_env$metacommunity,  by="terms", data = week8phyto_env)

phyto_perman12 <-adonis2(decostand(week12phyto,method="hellinger") ~dispersal*temp,  method = "bray", strata = week12phyto_env$metacommunity,  by="terms", data = week12phyto_env)


### PERMDISP

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

#Response: DistToCentroid
#          Df  Sum Sq  Mean Sq F value  Pr(>F)  
#dispersal  2 0.11356 0.056780  5.0124 0.01129 *
#Residuals 41 0.46444 0.011328     


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

# get post-hoc significance letters
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

phyletters_df

phyletters_letters_df <- phyletters_letters_df %>% 
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


# BACTERIA ----------------------------------------------------------------

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

#bacteria_perman4 <-adonis2(week4bacteria ~dispersal+temp_mean,  method = "bray", strata = week4bacteria_env$metacommunity,  by="terms", data = week4bacteria_env)
#bacteria_perman12 <-adonis2(week12bacteria ~dispersal+temp_mean,  method = "bray", strata = week12bacteria_env$metacommunity,  by="terms", data = week12bacteria_env)

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

# -------------------------------------------------------------------
### Now combine all for manuscript Fig. 6


disp_pal <- c("#003f5c", "#ff6361", "#ffa600")

LEG_TXT  <- 13.5
LEG_TTL  <- 14.5
XTICK    <- 13.5
XTITLE   <- 14.5
TITLE_SZ <- 14
YLAB     <- "Distance to metacommunity centroid"

common_fill <- scale_fill_manual(values = disp_pal, name = "Dispersal", limits = c("none","low","high"))


title_theme <- theme(
  plot.title.position = "plot",
  plot.title = element_text(size = TITLE_SZ, hjust = 0.5, margin = margin(l = 8))
)

# Phyto size structure 
p_size <- ggplot(
  phyto_size_bd_all,
  aes(x = factor(week, levels = c(2,4,6,8,10,12)),
      y = DistToCentroid, fill = dispersal)
) +
  geom_boxplot(width = 0.80, outlier.shape = NA, outlier.stroke = 0)+
  common_fill +
  scale_y_continuous(limits = c(0, 0.65), breaks = c(0.2,0.6), expand = c(0, 0)) +
  scale_x_discrete(name = "", labels = c("2","4","6","8","10","12")) +
  labs(y="") +
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
    plot.margin     = margin(22, 6, 5, 10)
  )

# Phyto taxonomic structure 
p_phyto <- ggplot(
  phyto_bd_all,
  aes(x = factor(week, levels = c(4,8,12)),
      y = DistToCentroid, fill = dispersal)
) +
  geom_boxplot(width = 0.7, outlier.shape = NA, outlier.stroke = 0) +
  common_fill +
  scale_y_continuous(limits = c(0.15, 0.75), breaks = c(0.3, 0.6), expand = c(0, 0)) +
  scale_x_discrete(name = "", labels = c("4","8","12")) +
  labs( y = "") +
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

# Bacterial taxonomic structure 
p_bacteria <- ggplot(
  bacteria_bd_all,
  aes(x = factor(week, levels = c(4,8,12)),
      y = DistToCentroid, fill = dispersal)
) +
  geom_boxplot(width = 0.7, outlier.shape = NA, outlier.stroke = 0) +
  common_fill +
  scale_y_continuous(limits = c(0.2, 0.6), breaks = c(0.3,0.5), expand = c(0, 0)) +
  scale_x_discrete(name = "", labels = c("4","8","12")) +
  labs( y = "") +
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

# Zoop taxonomic structure 
p_zoop <- ggplot(
  zoop_bd_all,
  aes(x = factor(week, levels = c(4,8,12)),
      y = DistToCentroid, fill = dispersal)
) +
  geom_boxplot(width = 0.7, outlier.shape = NA, outlier.stroke = 0) +
  common_fill +
  scale_y_continuous(limits = c(0, 0.8), breaks = c(0.2,0.6), expand = c(0, 0)) +
  scale_x_discrete(name = "", labels = c("4","8","12")) +
  labs( y = "") +
  theme_classic(base_size = 13) +
  title_theme +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = LEG_TXT),
    legend.title    = element_text(size = LEG_TTL),
    panel.border    = element_rect(colour = "black", fill = NA),
    axis.text.x     = element_text(size = XTICK),
    axis.text.y     = element_text(size = 12),
    axis.title.x    = element_text(size = XTITLE,vjust=1),
    plot.margin     = margin(6, 10, 10, 10)   
  )

# Assemble with custom widths since diff # of sampling times per panel
top_row    <- (p_size | p_bacteria) + plot_layout(widths = c(2.5, 0.9))
bottom_row <- (p_phyto | p_zoop)   + plot_layout(widths = c(1, 1))

final_permdisp_plot <-
  (top_row / bottom_row) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.margin = margin(15,15,15,15)  
  ) &
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


final_permdisp_plot +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(
    plot.tag.position = c(0.01, 1),  
    plot.tag = element_text(size = 14, face = "bold")
  )

#ggsave(final_permdisp_plot,"permdisp_final.pdf",width=8.4,height=7.8)