# Supplementary Figures

####  NUTRIENT CONCENTRATIONS OVER TIME
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


# Phytoplankton biomass in each individual mesocosm over time, separated by wattage treatment and dispersal 
biomass_disp <- ggplot(exp_data, aes(x = week_number,  y = ln_fluoro_ug_C_L,group=interaction(tank,dispersal)  ,color = dispersal,fill=dispersal)) + geom_point() +geom_line()+facet_grid(dispersal~wattage)+
  scale_color_manual(
    values = c("none" = "#003f5c", "low" = "#ff6361", "high" = "#ffa600"),
    labels = c("none" = "None", "low" = "Low", "high" = "High")
  )+
  scale_fill_manual(
    values = c("none" = "#003f5c", "low" = "#ff6361", "high" = "#ffa600"),
    labels = c("none" = "None", "low" = "Low", "high" = "High")
  )+
  theme_classic() +
  labs(
    x = "Week",
    y = expression("Phyto biomass (µg C L"^-1 * ")")
  )+ scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5)   # 1 decimal place
  ) +
  theme(legend.position= "none",
        strip.background = element_rect(fill = "white", color = "black"),
        panel.background = element_blank(),
        panel.border     = element_rect(color = "black", fill = NA),
        plot.background  = element_blank(),
        panel.grid       = element_blank(),
        strip.text     =element_text(size = 14, family = "Helvetica"),
        axis.title.x     = element_text(size = 14, family = "Helvetica",
                                        margin = margin(t = 11)),
        axis.title.y     = element_text(size = 14, family = "Helvetica",
                                        margin = margin(r = 11)),
        axis.text.x      = element_text(size = 13, family = "Helvetica"),
        axis.text.y      = element_text(size = 12, family = "Helvetica"),
        legend.text      = element_text(size = 13, family = "Helvetica"),
        legend.title     = element_text(size = 14, family = "Helvetica"),
        legend.key.width = unit(1, "cm")
  )

biomass_disp
ggsave(biomass_disp,file="biom_grid.pdf",width=8.5,height=6)



############## Visualizing variance patterns over time

# Summarize variance within each week × wattage × dispersal
var_summary_gpp <- exp_data %>%
  group_by(week_number, dispersal) %>%
  summarise(var_gpp = sd(GPP_umol_L_h, na.rm = TRUE),
            .groups = "drop")

gpp_variance_plot <- ggplot(var_summary_gpp, aes(x = dispersal, y = var_gpp, fill = dispersal)) +
  geom_boxplot()+
  scale_fill_manual(
    values = c("none" = "#003f5c", "low" = "#ff6361", "high" = "#ffa600"),
    labels = c("none" = "None", "low" = "Low", "high" = "High")
  ) +
  theme_classic() +
  labs(
    x = "",
    y = expression(STDEV~GPP~"(µmol L h"^"-1"*")")
  ) +
  scale_y_continuous(
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
                                    margin = margin(r = 11)),
    axis.text.x      = element_text(size = 13, family = "Helvetica"),
    axis.text.y      = element_text(size = 12, family = "Helvetica"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 13, family = "Helvetica"),
    legend.title     = element_text(size = 14, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(1, "cm")
  )

var_summary_biom <- exp_data %>%
  group_by(week_number, dispersal) %>%
  summarise(var_fluoro = sd(total_fluoro_ug_C_L, na.rm = TRUE),
            .groups = "drop")

biomass_variance_barplot <- ggplot(var_summary_biom, aes(x = dispersal, 
                                                         y = var_fluoro, 
                                                         fill = dispersal)) +
  geom_boxplot() +
  scale_fill_manual(
    values = c("none" = "#003f5c", "low" = "#ff6361", "high" = "#ffa600"),
    labels = c("none" = "None", "low" = "Low", "high" = "High")
  )+
  theme_classic() +
  labs(
    x = "Week",
    y = expression("STDEV Phyto biomass (µg C L"^-1 * ")")
  )+ scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::number_format(accuracy = 0.1)   # 1 decimal place
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
                                    margin = margin(r = 11),hjust=1),
    axis.text.x      = element_text(size = 13, family = "Helvetica"),
    axis.text.y      = element_text(size = 12, family = "Helvetica"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 13, family = "Helvetica"),
    legend.title     = element_text(size = 14, family = "Helvetica"),
    legend.margin    = margin(t = -6, b = 0, unit = "pt"),
    legend.key.width = unit(1, "cm")
  )

legend_var <- ggpubr::get_legend(
  gpp_variance_barplot + theme(legend.position = "bottom")
)

biom_var_nolegend <- biomass_variance_barplot + theme(legend.position="none",
                                                      plot.margin  = margin(t = 0, r = 5, b = 0, l = 5))

gpp_var_nolegend <- gpp_variance_barplot + theme(legend.position="none",
                                                 plot.margin  = margin(t = 5, r = 5, b = 0, l = 5))
phyto_alpha_div_plot<- ggplot(exp_data,
                              aes(x      = as.numeric(as.character(week_number)),
                                  y      = phyto_richness,
                                  colour = dispersal,
                                  group  = dispersal)) +
  stat_summary(fun.data = mean_se,
               geom     = "errorbar",
               position = pos,
               width    = 0.35,
               size     = 0.6,
               na.rm    = TRUE) +
  stat_summary(fun = mean,
               geom = "line",
               position = pos,
               size = 1,
               na.rm = TRUE) +
  
  ## mean points
  stat_summary(fun = mean,
               geom = "point",
               position = pos,
               size = 3,
               na.rm = TRUE) +
  
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
  scale_colour_manual(values = c("none" = "#003f5c",
                                 "low"  = "#ff6361",
                                 "high" = "#ffa600"),
                      name = "Dispersal") +
  labs(x = "Week number",
       y = "Phytoplankton alpha diversity") + 
  theme(legend.position = "right",
        legend.key.size     = unit(1, "cm"),
        panel.background    = element_blank(),
        panel.border        = element_rect(color = "black", fill = NA),
        plot.background     = element_blank(),
        panel.grid          = element_blank(),
        strip.text          = element_text(size = 16, face = "italic"),
        axis.title.x        = element_text(size = 16, family = "Helvetica"),
        axis.title.y        = element_text(size = 16, family = "Helvetica"),
        axis.text.x         = element_text(size = 14, family = "Helvetica"),
        axis.text.y         = element_text(size = 10, family = "Helvetica"),
        legend.text         = element_text(size = 10, family = "Helvetica"),
        legend.title        = element_text(size = 13, family = "Helvetica")
  ) 

zoop_alpha_div_plot<- ggplot(exp_data,
                             aes(x      = as.numeric(as.character(week_number)),
                                 y      = zoop_richness,
                                 colour = dispersal,
                                 group  = dispersal)) +
  
  ## error bars (mean ± SE)
  stat_summary(fun.data = mean_se,
               geom     = "errorbar",
               position = pos,
               width    = 0.35,
               size     = 0.6,
               na.rm    = TRUE) +
  
  ## means connected by lines
  stat_summary(fun = mean,
               geom = "line",
               position = pos,
               size = 1,
               na.rm = TRUE) +
  
  ## mean points
  stat_summary(fun = mean,
               geom = "point",
               position = pos,
               size = 3,
               na.rm = TRUE) +
  
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
  scale_colour_manual(values = c("none" = "#003f5c",
                                 "low"  = "#ff6361",
                                 "high" = "#ffa600"),
                      name = "Dispersal") +
  labs(x = "Week number",
       y = "Zooplankton alpha diversity") + 
  theme(legend.position = "right",
        legend.key.size     = unit(1, "cm"),
        panel.background    = element_blank(),
        panel.border        = element_rect(color = "black", fill = NA),
        plot.background     = element_blank(),
        panel.grid          = element_blank(),
        strip.text          = element_text(size = 16, face = "italic"),
        axis.title.x        = element_text(size = 16, family = "Helvetica"),
        axis.title.y        = element_text(size = 16, family = "Helvetica"),
        axis.text.x         = element_text(size = 14, family = "Helvetica"),
        axis.text.y         = element_text(size = 10, family = "Helvetica"),
        legend.text         = element_text(size = 10, family = "Helvetica"),
        legend.title        = element_text(size = 13, family = "Helvetica")
  ) 
