# Cleaning productivity data

install.packages("LakeMetabolizer")
library(LakeMetabolizer)
library(hms)
library(purrr)
prod<-read.csv("TEAMM_Productivity_trimmed.csv")

prod <- prod %>% mutate(tank = gsub("T", "", Tank))

prod <- prod %>%
  filter(tank != "Air", !is.na(tank)) %>%
  mutate(
    sample_date = as_date(
      parse_date_time(str_trim(sample_date), orders = c("mdY", "mdy"))
    )
  )

print(unique(prod$sample_date))
# Should give: 2018-06-27 2018-07-11 2018-07-25 2018-08-08 2018-08-22 2018-09-05

wind_data <- data.frame(
  sample_date_parsed = as_date(c("2018-06-27","2018-07-11","2018-07-25",
                                 "2018-08-08","2018-08-22","2018-09-05")),
  wind_kmh_airport   = c(14, 9, 10, 10, 11, 5),
  stringsAsFactors   = FALSE
) %>%
  mutate(
    wind_ms_site = (wind_kmh_airport / 3.6) * SHELTER_BEST,
    photoperiod  = sapply(sample_date_parsed, photoperiod_h)
  )

cat("\nwind_data:\n")
print(wind_data)

# ============================================================
# Now run the rest of the workflow unchanged from Step 1 onward
# replacing every instance of sample_date with sample_date_parsed
# ============================================================

prod_corrected <- prod %>%
  mutate(
    # Extract just the time part after the space - much more robust
    # Timestamp format is "M/D/YY H:MM" so split on space and take second element
    time_str = str_trim(str_extract(Timestamp, "\\d+:\\d+(:\\d+)?")),
    time     = as_hms(paste0(
      # pad to HH:MM:SS format
      str_pad(str_extract(time_str, "^\\d+"), 2, pad = "0"), ":",
      str_extract(time_str, "(?<=:)\\d+"), ":00"
    ))
  ) %>%
  left_join(wind_data, by = "sample_date") %>%
  mutate(
    do_sat_actual       = do_sat_mgL(temp),
    deficit_actual      = dissolved_O2_mg_L - do_sat_actual,
    k_o2_actual         = k_o2_md(wind_ms_site, temp),
    reair_rate_mg_L_day = k_o2_actual * deficit_actual / tank_depth_m
  )

# Check parse worked
cat("NA times:", sum(is.na(prod_corrected$time)), "\n")
cat("Sample times:\n")
print(head(prod_corrected %>% select(Timestamp, sample_date_parsed, time_of_day, time), 12))

cat("\nNAs in wind join:", sum(is.na(prod_corrected$wind_ms_site)), "\n")
cat("Rows:", nrow(prod_corrected), "\n")

# Window means
window_mean <- prod_corrected %>%
  filter(time_of_day %in% c("dawn1", "dusk", "dawn2")) %>%
  group_by(sample_date_parsed, tank, time_of_day) %>%
  summarise(
    dissolved_O2_mg_L   = mean(dissolved_O2_mg_L,     na.rm = TRUE),
    temp                = mean(temp,                   na.rm = TRUE),
    do_sat_actual       = mean(do_sat_actual,          na.rm = TRUE),
    deficit_actual      = mean(deficit_actual,         na.rm = TRUE),
    reair_rate_mg_L_day = mean(reair_rate_mg_L_day,   na.rm = TRUE),
    time                = mean(time,                   na.rm = TRUE),
    wind_ms_site        = first(wind_ms_site),
    photoperiod         = first(photoperiod),
    .groups = "drop"
  )

cat("\nwindow_mean rows:", nrow(window_mean), "\n")


prod$Cstar <- o2.at.sat.base(prod$temp,            
                        baro = 1013.25) 
str(prod)
T_ref       <- 23                                  
Cstar_ref   <- o2.at.sat.base(T_ref, baro = 1013.25)   


prod <- prod %>%
  mutate(Timestamp   = str_trim(Timestamp),      
         datetime    = parse_date_time(
           Timestamp,
           orders = c("mdY HMS",   
                      "mdY HM",    
                      "dmy HMS",   
                      "dmy HM"),
           tz = "UTC")             
  )

prod <- prod %>%                         
         mutate(DO_temp_adj    = dissolved_O2_mg_L * Cstar_ref / Cstar) 

prod <- prod %>%
  mutate(
    date = as_date(datetime),       
    time = as_hms(datetime)          
   
  )

#  Collapse each 30‑min window to a single mean value 

window_mean <- prod %>%                       
  group_by(sample_date, tank, time_of_day) %>%
  dplyr::summarise(
    DO_temp_adj = mean(DO_temp_adj, na.rm = TRUE),
    time        = mean(time,         na.rm = TRUE), 
    .groups     = "keep"
  )

cycle <- window_mean %>%
  pivot_wider(
    id_cols     = c(sample_date, tank),
    names_from  = time_of_day,
    values_from = c(DO_temp_adj, time),
    names_sep   = "_"
  )

sec_day <- 24 * 3600 #convert day to seconds

cycle <- cycle %>%
  mutate(
    # keep only clock time portion
    time_dawn1 = as_hms(time_dawn1),
    time_dusk  = as_hms(time_dusk),
    time_dawn2 = as_hms(time_dawn2),
    
    # daylight duration 
    hours_day = as.numeric((time_dusk - time_dawn1) , units = "secs") / 3600,
    
    # nighttime duration: wrap negative gaps across midnight
    sec_night = as.numeric((time_dawn2 - time_dusk) , units = "secs"),
    hours_night = if_else(sec_night < 0,
                          (sec_night + sec_day) / 3600,
                          sec_night          / 3600),
  
    NPP = (DO_temp_adj_dusk  - DO_temp_adj_dawn1) / hours_day,
    ER  = (DO_temp_adj_dawn2 - DO_temp_adj_dusk)  / hours_night
  ) %>%
  select(sample_date, tank, NPP, ER) 

fluxes_clean <- cycle %>%mutate(
    NPP_umol = NPP * 31.25,  # convert mg O2 L/h to umol
    ER_umol  = ER  * 31.25  ,
    GPP_umol = NPP_umol-ER_umol
  )
