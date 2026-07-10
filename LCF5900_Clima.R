install.packages(c(
  "rio",           # Universal file import/export
  "tidyverse",     # dplyr, ggplot2, tidyr, etc.
  "lubridate",     # Date manipulation
  "broom",         # Tidy model outputs
  "zoo",           # Rolling means
  "trend",         # Mann-Kendall and Sen's slope
  "scales",        # Axis formatting
  "patchwork"      # Combine ggplots))
))

  library(rio)
  library(tidyverse)
  library(lubridate)
  library(broom)
  library(zoo)
  library(trend)
  library(scales)
  library(patchwork)
  
  # Import the full dataset from Excel
  # The rio package auto-detects the format from the file extension
  climate_data <- import("DadosClima_Piracicaba.xlsx")
  
  # Create a proper Date column and set factor levels
  climate_data <- climate_data %>%
    mutate(
      # Create Date object
      Date = make_date(Ano, Mes, Dia),
      
      # Define Trimestre as ordered factor (rolling 3-month windows)
      Trimestre = factor(Trimestre,
                         levels = c("JFM", "FMA", "MAM", "AMJ", "MJJ",
                                    "JJA", "JAS", "ASO", "SON", "OND", "NDJ", "DJF")),
      
      # ENSO classification
      ClassNino = factor(ClassNino,
                         levels = c("StrgNina", "Interm", "StrgNino"),
                         labels = c("Strong La Ni\u00f1a", "Intermediate", "Strong El Ni\u00f1o")),
      
      # Seasons (Southern Hemisphere)
      Estacao = factor(Estacao,
                       levels = c("Ver\u00e3o", "Outono", "Inverno", "Primavera"),
                       labels = c("Summer", "Autumn", "Winter", "Spring")),
      
      # Decade for grouping analyses
      Decade = floor(Ano / 10) * 10,
      
      # Diurnal Temperature Range
      DTR = TMAX - TMIN
    )
  
  # Quick summary
  cat("Dataset spans from", min(climate_data$Date, na.rm = TRUE),
      "to", max(climate_data$Date, na.rm = TRUE), "\n")
  cat("Total observations:", nrow(climate_data), "\n")
  
  # Function: Mann-Kendall trend test with Sen's slope
  trend_summary <- function(x, y) {
    # Remove NAs
    valid <- complete.cases(x, y)
    x <- x[valid]; y <- y[valid]
    
    # Mann-Kendall test
    mk <- mk.test(y)
    
    # Sen's slope
    ss <- sens.slope(y)
    
    tibble(
      Sen_Slope = ss$estimates,
      MK_Tau = mk$estimates,
      p_value = mk$p.value,
      Significant = ifelse(mk$p.value < 0.05, "Yes", "No")
    )
  }
  
  # Function: Calculate moving average
  moving_avg <- function(x, n = 10) {
    zoo::rollmean(x, k = n, fill = NA, align = "center")
  }
  
  # Calculate annual means for Tmin, Tmed, Tmax
  annual_temp <- climate_data %>%
    filter(Ano >= 1917, !is.na(TMED), !is.na(TMIN), !is.na(TMAX)) %>%
    group_by(Ano) %>%
    summarise(
      Mean_TMIN = mean(TMIN, na.rm = TRUE),
      Mean_TMED = mean(TMED, na.rm = TRUE),
      Mean_TMAX = mean(TMAX, na.rm = TRUE),
      Mean_DTR  = mean(TMAX - TMIN, na.rm = TRUE),
      .groups = "drop"
    )
  
  # View the first rows
  head(annual_temp)
  
  # Linear models
  lm_tmin <- lm(Mean_TMIN ~ Ano, data = annual_temp)
  lm_tmed <- lm(Mean_TMED ~ Ano, data = annual_temp)
  lm_tmax <- lm(Mean_TMAX ~ Ano, data = annual_temp)
  
  # Extract slopes (warming rate in °C per year)
  slopes <- tibble(
    Variable = c("TMIN", "TMED", "TMAX"),
    Slope_per_year = c(coef(lm_tmin)[2], coef(lm_tmed)[2], coef(lm_tmax)[2]),
    Slope_per_decade = Slope_per_year * 10,
    R_squared = c(summary(lm_tmin)$r.squared,
                  summary(lm_tmed)$r.squared,
                  summary(lm_tmax)$r.squared),
    p_value = c(summary(lm_tmin)$coefficients[2,4],
                summary(lm_tmed)$coefficients[2,4],
                summary(lm_tmax)$coefficients[2,4])
  )
  
  print(slopes)
  
  # Mann-Kendall non-parametric trend test
  cat("\n--- Mann-Kendall Trend Tests ---\n")
  cat("TMIN: "); print(trend_summary(annual_temp$Ano, annual_temp$Mean_TMIN))
  cat("TMED: "); print(trend_summary(annual_temp$Ano, annual_temp$Mean_TMED))
  cat("TMAX: "); print(trend_summary(annual_temp$Ano, annual_temp$Mean_TMAX))
  
  # Reshape for plotting
  temp_long <- annual_temp %>%
    pivot_longer(cols = c(Mean_TMIN, Mean_TMED, Mean_TMAX),
                 names_to = "Variable", values_to = "Temperature") %>%
    mutate(Variable = recode(Variable,
                             "Mean_TMIN" = "Tmin", "Mean_TMED" = "Tmean", "Mean_TMAX" = "Tmax"))

  # Plot with linear trends
  p1 <- ggplot(temp_long, aes(x = Ano, y = Temperature, color = Variable)) +
    geom_line(alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
    scale_color_manual(values = c("Tmin" = "#3498db", "Tmean" = "#2ecc71", "Tmax" = "#e74c3c")) +
    labs(
      title = "Long-term Temperature Trends in Piracicaba (1917\u2013Present)",
      subtitle = "Linear trends fitted to annual means",
      x = "Year", y = "Temperature (\u00b0C)", color = "Variable"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  print(p1)
  ggsave("Q1_temperature_trends.png", p1, width = 10, height = 6, dpi = 150)
  
  # DTR trend analysis
  lm_dtr <- lm(Mean_DTR ~ Ano, data = annual_temp)
  cat("DTR trend (°C/decade):", coef(lm_dtr)[2] * 10, "\n")
  cat("DTR Mann-Kendall:\n")
  print(trend_summary(annual_temp$Ano, annual_temp$Mean_DTR))
  
  # Plot DTR
  p2 <- ggplot(annual_temp, aes(x = Ano, y = Mean_DTR)) +
    geom_line(color = "grey50") +
    geom_smooth(method = "lm", color = "#e74c3c", fill = "#fadbd8") +
    geom_line(aes(y = moving_avg(Mean_DTR, 10)), color = "#2c3e50", linewidth = 1) +
    labs(
      title = "Diurnal Temperature Range (Tmax \u2212 Tmin)",
      subtitle = "Annual mean with 10-year moving average (black) and linear trend (red)",
      x = "Year", y = "DTR (\u00b0C)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p2)
  ggsave("Q1_DTR_trend.png", p2, width = 10, height = 5, dpi = 150)
  
  # Count days exceeding 35°C per year
  extreme_heat <- climate_data %>%
    filter(Ano >= 1917, !is.na(TMAX)) %>%
    group_by(Ano) %>%
    summarise(
      Days_Above_35 = sum(TMAX > 35, na.rm = TRUE),
      Days_Above_33 = sum(TMAX > 33, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Trend test
  cat("Trend in days > 35°C:\n")
  print(trend_summary(extreme_heat$Ano, extreme_heat$Days_Above_35))
  
  # Plot
  p3 <- ggplot(extreme_heat, aes(x = Ano, y = Days_Above_35)) +
    geom_col(fill = "#c0392b", alpha = 0.7) +
    geom_smooth(method = "loess", span = 0.3, color = "black", se = FALSE) +
    labs(
      title = "Annual Frequency of Extreme Heat Days (Tmax > 35\u00b0C)",
      x = "Year", y = "Number of Days"
    ) +
    theme_minimal(base_size = 12)
  
  print(p3)
  ggsave("Q1_extreme_heat_days.png", p3, width = 10, height = 5, dpi = 150)
  
  # Detect structural breaks in the temperature series
  # Using the 'strucchange' package
  # install.packages("strucchange")
  library(strucchange)
  
  # Test for breakpoints in TMIN annual series
  bp_tmin <- breakpoints(Mean_TMIN ~ Ano, data = annual_temp, h = 15)
  summary(bp_tmin)
  
  # Plot with breakpoints
  plot(bp_tmin)
  cat("Breakpoint years for TMIN:", annual_temp$Ano[bp_tmin$breakpoints], "\n")
  
  # Annual total rainfall
  annual_rain <- climate_data %>%
    filter(Ano >= 1902, !is.na(Chuva)) %>%
    group_by(Ano) %>%
    summarise(
      Total_Rain = sum(Chuva, na.rm = TRUE),
      Rainy_Days = sum(Chuva > 0.1, na.rm = TRUE),
      Mean_Intensity = Total_Rain / Rainy_Days,
      .groups = "drop"
    )
  
  # Mann-Kendall trend test
  cat("Annual rainfall trend:\n")
  print(trend_summary(annual_rain$Ano, annual_rain$Total_Rain))
  
  # Plot
  p4 <- ggplot(annual_rain, aes(x = Ano, y = Total_Rain)) +
    geom_line(color = "#2980b9", alpha = 0.6) +
    geom_smooth(method = "lm", color = "#e74c3c", fill = "#fadbd8") +
    geom_line(aes(y = moving_avg(Total_Rain, 10)), color = "#2c3e50", linewidth = 1) +
    geom_hline(yintercept = mean(annual_rain$Total_Rain), linetype = "dashed", color = "grey50") +
    labs(
      title = "Total Annual Rainfall in Piracicaba (1902\u2013Present)",
      subtitle = "10-year moving average (black), linear trend (red), long-term mean (dashed)",
      x = "Year", y = "Annual Rainfall (mm)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p4)
  ggsave("Q2_annual_rainfall.png", p4, width = 10, height = 5, dpi = 150)
  
  # Maximum dry spell per year
  dry_spells <- climate_data %>%
    filter(Ano >= 1902) %>%
    group_by(Ano) %>%
    summarise(
      Max_Estiagem = max(Estiagem, na.rm = TRUE),
      Mean_Estiagem = mean(Estiagem[Estiagem > 0], na.rm = TRUE),
      N_Long_Spells = sum(Estiagem > 20, na.rm = TRUE),  # Days in spells > 20 days
      .groups = "drop"
    )
  
  # Trend tests
  cat("Max dry spell trend:\n")
  print(trend_summary(dry_spells$Ano, dry_spells$Max_Estiagem))
  
  # Plot maximum dry spell length over time
  p5 <- ggplot(dry_spells, aes(x = Ano, y = Max_Estiagem)) +
    geom_line(color = "#d35400", alpha = 0.6) +
    geom_smooth(method = "lm", color = "#c0392b") +
    geom_line(aes(y = moving_avg(Max_Estiagem, 10)), color = "#2c3e50", linewidth = 1) +
    labs(
      title = "Maximum Annual Dry Spell Length",
      subtitle = "Longest consecutive period without rain each year",
      x = "Year", y = "Days"
    ) +
    theme_minimal(base_size = 12)
  
  print(p5)
  ggsave("Q2_dry_spells.png", p5, width = 10, height = 5, dpi = 150)
  
  # Decadal comparison boxplot
  dry_spell_decade <- climate_data %>%
    filter(Ano >= 1902, Estiagem > 0) %>%
    mutate(Decade = factor(floor(Ano / 10) * 10))
  
  p5b <- ggplot(dry_spell_decade, aes(x = Decade, y = Estiagem)) +
    geom_boxplot(fill = "#f39c12", alpha = 0.5, outlier.alpha = 0.3) +
    labs(
      title = "Distribution of Dry Spell Lengths by Decade",
      x = "Decade", y = "Dry Spell Length (days)"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p5b)
  ggsave("Q2_dry_spells_decades.png", p5b, width = 10, height = 5, dpi = 150)
  
  extreme_rain <- climate_data %>%
    filter(Ano >= 1902, !is.na(Chuva)) %>%
    group_by(Ano) %>%
    summarise(
      Max_Daily = max(Chuva, na.rm = TRUE),
      P95 = quantile(Chuva[Chuva > 0.1], 0.95, na.rm = TRUE),
      P99 = quantile(Chuva[Chuva > 0.1], 0.99, na.rm = TRUE),
      Days_Above_50mm = sum(Chuva > 50, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Trend in annual maximum daily rainfall
  cat("Max daily rainfall trend:\n")
  print(trend_summary(extreme_rain$Ano, extreme_rain$Max_Daily))
  
  # Plot
  p6 <- ggplot(extreme_rain, aes(x = Ano, y = Max_Daily)) +
    geom_point(color = "#2980b9", alpha = 0.6) +
    geom_smooth(method = "lm", color = "#e74c3c") +
    labs(
      title = "Annual Maximum Daily Rainfall",
      x = "Year", y = "Maximum Daily Rainfall (mm)"
    ) +
    theme_minimal(base_size = 12)
  
  # Days above 50mm
  p6b <- ggplot(extreme_rain, aes(x = Ano, y = Days_Above_50mm)) +
    geom_col(fill = "#2980b9", alpha = 0.6) +
    geom_smooth(method = "loess", span = 0.3, color = "black", se = FALSE) +
    labs(
      title = "Annual Count of Heavy Rain Days (> 50 mm/day)",
      x = "Year", y = "Number of Days"
    ) +
    theme_minimal(base_size = 12)
  
  print(p6 / p6b)
  ggsave("Q2_extreme_rainfall.png", p6 / p6b, width = 10, height = 8, dpi = 150)
  
  # Define wet season onset as the first day after September 1
  # with accumulated rainfall > 20mm in a 5-day window
  wet_season_onset <- climate_data %>%
    filter(!is.na(Chuva), Mes %in% c(9, 10, 11, 12)) %>%
    arrange(Date) %>%
    group_by(Ano) %>%
    mutate(
      Rain_5day = zoo::rollsum(Chuva, k = 5, fill = NA, align = "right")
    ) %>%
    filter(Rain_5day >= 20) %>%
    slice_min(Date, n = 1) %>%
    ungroup() %>%
    mutate(
      # Day of year for the onset
      DOY_onset = yday(Date)
    ) %>%
    select(Ano, Date, DOY_onset)
  
  # Trend in onset date
  cat("Wet season onset trend:\n")
  print(trend_summary(wet_season_onset$Ano, wet_season_onset$DOY_onset))
  
  # Plot
  p8 <- ggplot(wet_season_onset, aes(x = Ano, y = DOY_onset)) +
    geom_point(color = "#27ae60", alpha = 0.6) +
    geom_smooth(method = "lm", color = "#c0392b") +
    scale_y_continuous(breaks = c(244, 274, 305, 335),
                       labels = c("Sep 1", "Oct 1", "Nov 1", "Dec 1")) +
    labs(
      title = "Wet Season Onset Date (1902\u2013Present)",
      subtitle = "First 5-day period with \u226520 mm after September 1",
      x = "Year", y = "Onset Date"
    ) +
    theme_minimal(base_size = 12)
  
  print(p8)
  ggsave("Q3_wet_season_onset.png", p8, width = 10, height = 5, dpi = 150)
  
  # Frost risk: days with Tmin < 3°C (severe frost) or Tmin < 5°C (light frost)
  frost_risk <- climate_data %>%
    filter(Ano >= 1917, !is.na(TMIN)) %>%
    group_by(Ano) %>%
    summarise(
      Frost_Severe = sum(TMIN < 3, na.rm = TRUE),
      Frost_Light  = sum(TMIN < 5, na.rm = TRUE),
      Last_Frost_DOY = ifelse(any(TMIN < 5),
                              max(yday(Date[TMIN < 5])), NA_real_),
      First_Frost_DOY = ifelse(any(TMIN < 5),
                               min(yday(Date[TMIN < 5])), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(Frost_Window = Last_Frost_DOY - First_Frost_DOY)
  
  # Trend in frost days
  cat("Light frost days trend:\n")
  print(trend_summary(frost_risk$Ano, frost_risk$Frost_Light))
  
  # Plot frost frequency
  p9 <- ggplot(frost_risk, aes(x = Ano, y = Frost_Light)) +
    geom_col(fill = "#85c1e9", alpha = 0.7) +
    geom_smooth(method = "loess", span = 0.3, color = "#2c3e50", se = FALSE) +
    labs(
      title = "Annual Frost Risk Days (Tmin < 5\u00b0C)",
      subtitle = "Is the frost-risk window shrinking?",
      x = "Year", y = "Number of Days"
    ) +
    theme_minimal(base_size = 12)
  
  print(p9)
  ggsave("Q3_frost_risk.png", p9, width = 10, height = 5, dpi = 150)
  
  # Monthly rainfall by Trimestre across decades
  trimestre_rain <- climate_data %>%
    filter(Ano >= 1902, !is.na(Chuva)) %>%
    mutate(Period = case_when(
      Ano < 1950 ~ "1902\u20131949",
      Ano < 1980 ~ "1950\u20131979",
      Ano < 2000 ~ "1980\u20131999",
      TRUE ~ "2000\u2013Present"
    )) %>%
    group_by(Period, Trimestre) %>%
    summarise(Mean_Rain = mean(Chuva, na.rm = TRUE), .groups = "drop")
  
  p10 <- ggplot(trimestre_rain, aes(x = Trimestre, y = Mean_Rain, fill = Period)) +
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Mean Daily Rainfall by Trimestre Across Periods",
      x = "Trimestre (Rolling 3-month Window)", y = "Mean Daily Rainfall (mm)"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p10)
  ggsave("Q3_trimestre_rainfall.png", p10, width = 10, height = 5, dpi = 150)
  
  # Seasonal rainfall by ENSO phase
  enso_rain <- climate_data %>%
    filter(Ano >= 1950, !is.na(ClassNino), !is.na(Chuva)) %>%
    group_by(Ano, Estacao, ClassNino) %>%
    summarise(Seasonal_Rain = sum(Chuva, na.rm = TRUE), .groups = "drop")
  
  # Boxplot comparison
  p11 <- ggplot(enso_rain, aes(x = ClassNino, y = Seasonal_Rain, fill = ClassNino)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.4) +
    facet_wrap(~ Estacao, scales = "free_y") +
    scale_fill_manual(values = c(
      "Strong La Ni\u00f1a" = "#3498db",
      "Intermediate" = "#95a5a6",
      "Strong El Ni\u00f1o" = "#e74c3c"
    )) +
    labs(
      title = "Seasonal Rainfall by ENSO Phase (1950\u2013Present)",
      x = "ENSO Classification", y = "Seasonal Rainfall (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  
  print(p11)
  ggsave("Q4_enso_rainfall.png", p11, width = 11, height = 7, dpi = 150)
  
  # Statistical test: Kruskal-Wallis for each season
  enso_rain %>%
    group_by(Estacao) %>%
    summarise(
      KW_stat = kruskal.test(Seasonal_Rain ~ ClassNino)$statistic,
      p_value = kruskal.test(Seasonal_Rain ~ ClassNino)$p.value,
      .groups = "drop"
    ) %>%
    print()
  
  # Calculate monthly temperature anomalies relative to 1961-1990 baseline
  baseline_temp <- climate_data %>%
    filter(Ano >= 1961, Ano <= 1990) %>%
    group_by(Mes) %>%
    summarise(Clim_TMED = mean(TMED, na.rm = TRUE), .groups = "drop")
  
  enso_temp <- climate_data %>%
    filter(Ano >= 1950, !is.na(ClassNino), !is.na(TMED)) %>%
    left_join(baseline_temp, by = "Mes") %>%
    mutate(Temp_Anomaly = TMED - Clim_TMED) %>%
    group_by(Ano, ClassNino) %>%
    summarise(Mean_Anomaly = mean(Temp_Anomaly, na.rm = TRUE), .groups = "drop")
  
  # Plot: Are anomalies diverging over time by ENSO phase?
  p12 <- ggplot(enso_temp, aes(x = Ano, y = Mean_Anomaly, color = ClassNino)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = c(
      "Strong La Ni\u00f1a" = "#3498db",
      "Intermediate" = "#95a5a6",
      "Strong El Ni\u00f1o" = "#e74c3c"
    )) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = "Temperature Anomalies by ENSO Phase Over Time",
      subtitle = "Is the ENSO signal amplifying as background climate warms?",
      x = "Year", y = "Temperature Anomaly (\u00b0C)", color = "ENSO Phase"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  print(p12)
  ggsave("Q4_enso_temp_signal.png", p12, width = 10, height = 6, dpi = 150)
  
  # Monthly correlations between ONI and local climate
  oni_correlations <- climate_data %>%
    filter(Ano >= 1950, !is.na(MeanONI)) %>%
    group_by(Mes) %>%
    summarise(
      Cor_Rain = cor(MeanONI, Chuva, use = "complete.obs"),
      Cor_Tmed = cor(MeanONI, TMED, use = "complete.obs"),
      Cor_Tmin = cor(MeanONI, TMIN, use = "complete.obs"),
      .groups = "drop"
    )
  
  # Reshape and plot
  oni_cor_long <- oni_correlations %>%
    pivot_longer(-Mes, names_to = "Variable", values_to = "Correlation")
  
  p13 <- ggplot(oni_cor_long, aes(x = factor(Mes), y = Correlation, fill = Variable)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 0) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Monthly Correlation Between ONI and Local Climate Variables",
      x = "Month", y = "Pearson Correlation"
    ) +
    theme_minimal(base_size = 12)
  
  print(p13)
  ggsave("Q4_oni_correlations.png", p13, width = 10, height = 5, dpi = 150)
  
  # TR20: Days with Tmin > 20°C
  tr20 <- climate_data %>%
    filter(Ano >= 1917, !is.na(TMIN)) %>%
    group_by(Ano) %>%
    summarise(TR20 = sum(TMIN > 20, na.rm = TRUE), .groups = "drop")
  
  # Trend
  cat("TR20 (Tropical Nights) trend:\n")
  print(trend_summary(tr20$Ano, tr20$TR20))
  
  p14 <- ggplot(tr20, aes(x = Ano, y = TR20)) +
    geom_line(color = "#e74c3c", alpha = 0.6) +
    geom_smooth(method = "lm", color = "#2c3e50") +
    labs(
      title = "Tropical Nights (TR20): Days with Tmin > 20\u00b0C",
      x = "Year", y = "Number of Nights"
    ) +
    theme_minimal(base_size = 12)
  
  print(p14)
  ggsave("Q5_tropical_nights.png", p14, width = 10, height = 5, dpi = 150)
  
  
  # CDD: Maximum consecutive dry days per year
  # The Estiagem column directly gives us this information
  cdd <- climate_data %>%
    filter(Ano >= 1902) %>%
    group_by(Ano) %>%
    summarise(CDD = max(Estiagem, na.rm = TRUE), .groups = "drop")
  
  cat("CDD trend:\n")
  print(trend_summary(cdd$Ano, cdd$CDD))
  
  p15 <- ggplot(cdd, aes(x = Ano, y = CDD)) +
    geom_line(color = "#d35400", alpha = 0.6) +
    geom_smooth(method = "lm", color = "#c0392b") +
    labs(
      title = "Consecutive Dry Days (CDD) Index",
      subtitle = "Maximum annual dry spell length from Estiagem variable",
      x = "Year", y = "CDD (days)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p15)
  ggsave("Q5_CDD.png", p15, width = 10, height = 5, dpi = 150)
  
  # Calculate baseline percentiles (1961-1990, wet days only: Chuva >= 1mm)
  wet_days_baseline <- climate_data %>%
    filter(Ano >= 1961, Ano <= 1990, Chuva >= 1)
  
  p95 <- quantile(wet_days_baseline$Chuva, 0.95, na.rm = TRUE)
  p99 <- quantile(wet_days_baseline$Chuva, 0.99, na.rm = TRUE)
  cat("95th percentile threshold:", p95, "mm\n")
  cat("99th percentile threshold:", p99, "mm\n")
  
  # Annual R95p and R99p
  r95p <- climate_data %>%
    filter(Ano >= 1902, !is.na(Chuva)) %>%
    group_by(Ano) %>%
    summarise(
      R95p = sum(Chuva[Chuva > p95], na.rm = TRUE),
      R99p = sum(Chuva[Chuva > p99], na.rm = TRUE),
      SDII = sum(Chuva[Chuva >= 1], na.rm = TRUE) / sum(Chuva >= 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Trends
  cat("\nR95p trend:\n"); print(trend_summary(r95p$Ano, r95p$R95p))
  cat("R99p trend:\n"); print(trend_summary(r95p$Ano, r95p$R99p))
  cat("SDII trend:\n"); print(trend_summary(r95p$Ano, r95p$SDII))
  
  # Plot R95p
  p16 <- ggplot(r95p, aes(x = Ano, y = R95p)) +
    geom_col(fill = "#2980b9", alpha = 0.6) +
    geom_smooth(method = "loess", span = 0.3, color = "#c0392b", se = FALSE) +
    labs(
      title = "R95p: Total Rainfall from Very Wet Days (> 95th percentile)",
      x = "Year", y = "R95p (mm)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p16)
  ggsave("Q5_R95p.png", p16, width = 10, height = 5, dpi = 150)
  
  # Calculate 90th percentile thresholds from baseline (1961-1990)
  tmax_baseline <- climate_data %>%
    filter(Ano >= 1961, Ano <= 1990, !is.na(TMAX))
  tmin_baseline <- climate_data %>%
    filter(Ano >= 1961, Ano <= 1990, !is.na(TMIN))
  
  tx90_threshold <- quantile(tmax_baseline$TMAX, 0.90, na.rm = TRUE)
  tn90_threshold <- quantile(tmin_baseline$TMIN, 0.90, na.rm = TRUE)
  
  # Annual percentage of warm days/nights
  warm_indices <- climate_data %>%
    filter(Ano >= 1917, !is.na(TMAX), !is.na(TMIN)) %>%
    group_by(Ano) %>%
    summarise(
      TX90p = sum(TMAX > tx90_threshold) / n() * 100,
      TN90p = sum(TMIN > tn90_threshold) / n() * 100,
      .groups = "drop"
    )
  
  # Plot
  warm_long <- warm_indices %>%
    pivot_longer(-Ano, names_to = "Index", values_to = "Percent")
  
  p17 <- ggplot(warm_long, aes(x = Ano, y = Percent, color = Index)) +
    geom_line(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("TX90p" = "#e74c3c", "TN90p" = "#8e44ad")) +
    labs(
      title = "Warm Days (TX90p) and Warm Nights (TN90p)",
      subtitle = "Percentage of days exceeding 90th percentile; dashed line = expected 10%",
      x = "Year", y = "Percentage (%)", color = "Index"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  print(p17)
  ggsave("Q5_warm_days_nights.png", p17, width = 10, height = 5, dpi = 150)
  
  # Annual mean global radiation
  annual_rad <- climate_data %>%
    filter(Ano >= 2000, !is.na(GlobRad)) %>%
    group_by(Ano) %>%
    summarise(
      Mean_GlobRad = mean(GlobRad, na.rm = TRUE),
      Median_GlobRad = median(GlobRad, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Trend test
  cat("Global radiation trend (2000-present):\n")
  lm_rad <- lm(Mean_GlobRad ~ Ano, data = annual_rad)
  summary(lm_rad)
  cat("Trend (MJ/m²/day per decade):", coef(lm_rad)[2] * 10, "\n")
  
  # Plot
  p18 <- ggplot(annual_rad, aes(x = Ano, y = Mean_GlobRad)) +
    geom_line(color = "#f39c12", linewidth = 1) +
    geom_point(color = "#f39c12", size = 2) +
    geom_smooth(method = "lm", color = "#c0392b", fill = "#fadbd8") +
    labs(
      title = "Annual Mean Global Radiation (2000\u2013Present)",
      subtitle = "Solar brightening or dimming?",
      x = "Year", y = "Mean Global Radiation (MJ/m\u00b2/day)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p18)
  ggsave("Q6_radiation_trend.png", p18, width = 10, height = 5, dpi = 150)
  
  # Seasonal radiation trends
  seasonal_rad <- climate_data %>%
    filter(Ano >= 2000, !is.na(GlobRad)) %>%
    group_by(Ano, Estacao) %>%
    summarise(Mean_GlobRad = mean(GlobRad, na.rm = TRUE), .groups = "drop")
  
  p18b <- ggplot(seasonal_rad, aes(x = Ano, y = Mean_GlobRad, color = Estacao)) +
    geom_line() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = "Seasonal Global Radiation Trends",
      x = "Year", y = "Mean Global Radiation (MJ/m\u00b2/day)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p18b)
  ggsave("Q6_seasonal_radiation.png", p18b, width = 10, height = 5, dpi = 150)
  
  # Multiple regression: Does radiation explain temperature beyond the trend?
  rad_data <- climate_data %>%
    filter(Ano >= 2000, !is.na(GlobRad), !is.na(TMED))
  
  # Model 1: Temperature ~ Year only (thermal trend)
  mod1 <- lm(TMED ~ Ano, data = rad_data)
  
  # Model 2: Temperature ~ Year + Radiation
  mod2 <- lm(TMED ~ Ano + GlobRad, data = rad_data)
  
  # Model 3: Temperature ~ Year + Radiation + Month (seasonal control)
  mod3 <- lm(TMED ~ Ano + GlobRad + factor(Mes), data = rad_data)
  
  # Compare models
  cat("Model 1 (Year only) R²:", summary(mod1)$r.squared, "\n")
  cat("Model 2 (Year + Rad) R²:", summary(mod2)$r.squared, "\n")
  cat("Model 3 (Year + Rad + Month) R²:", summary(mod3)$r.squared, "\n")
  
  # ANOVA comparison
  anova(mod1, mod2, mod3)
  
  # Does cloud cover (implied by lower radiation) vary with ENSO?
  rad_enso <- climate_data %>%
    filter(Ano >= 2000, !is.na(GlobRad), !is.na(ClassNino)) %>%
    group_by(Estacao, ClassNino) %>%
    summarise(Mean_Rad = mean(GlobRad, na.rm = TRUE), .groups = "drop")
  
  p19 <- ggplot(rad_enso, aes(x = ClassNino, y = Mean_Rad, fill = ClassNino)) +
    geom_col(alpha = 0.7) +
    facet_wrap(~ Estacao) +
    scale_fill_manual(values = c(
      "Strong La Ni\u00f1a" = "#3498db",
      "Intermediate" = "#95a5a6",
      "Strong El Ni\u00f1o" = "#e74c3c"
    )) +
    labs(
      title = "Mean Global Radiation by ENSO Phase and Season",
      subtitle = "Lower radiation implies more cloud cover",
      x = "ENSO Phase", y = "Mean Radiation (MJ/m\u00b2/day)"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  
  print(p19)
  ggsave("Q6_radiation_enso.png", p19, width = 10, height = 6, dpi = 150)
  
  # Define compound event thresholds
  # Drought: Estiagem > 15 consecutive days
  # Heat: Tmax > 33°C (above the 90th percentile for the region)
  compound <- climate_data %>%
    filter(Ano >= 1917, !is.na(TMAX)) %>%
    mutate(
      Is_Drought = Estiagem > 15,
      Is_Heat = TMAX > 33,
      Is_Compound = Is_Drought & Is_Heat
    )
  
  # Annual count of compound days
  compound_annual <- compound %>%
    group_by(Ano) %>%
    summarise(
      Drought_Days = sum(Is_Drought, na.rm = TRUE),
      Heat_Days = sum(Is_Heat, na.rm = TRUE),
      Compound_Days = sum(Is_Compound, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Trend
  cat("Compound event days trend:\n")
  print(trend_summary(compound_annual$Ano, compound_annual$Compound_Days))
  
  # Plot
  p20 <- ggplot(compound_annual, aes(x = Ano, y = Compound_Days)) +
    geom_col(fill = "#8e44ad", alpha = 0.7) +
    geom_smooth(method = "loess", span = 0.3, color = "black", se = FALSE) +
    labs(
      title = "Compound Heat-Drought Days per Year",
      subtitle = "Days with Estiagem > 15 AND Tmax > 33\u00b0C",
      x = "Year", y = "Number of Compound Days"
    ) +
    theme_minimal(base_size = 12)
  
  print(p20)
  ggsave("Q7_compound_events.png", p20, width = 10, height = 5, dpi = 150)
  
  # Focus on dry season (Inverno/Winter) under different ENSO phases
  dry_season_compound <- climate_data %>%
    filter(Ano >= 1950, Estacao == "Winter", !is.na(ClassNino)) %>%
    group_by(Ano, ClassNino) %>%
    summarise(
      Mean_TMAX = mean(TMAX, na.rm = TRUE),
      Total_Rain = sum(Chuva, na.rm = TRUE),
      Max_Estiagem = max(Estiagem, na.rm = TRUE),
      Compound_Days = sum(Estiagem > 15 & TMAX > 33, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Scatter: Temperature vs Rainfall colored by ENSO
  p21 <- ggplot(dry_season_compound, aes(x = Total_Rain, y = Mean_TMAX, color = ClassNino)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = c(
      "Strong La Ni\u00f1a" = "#3498db",
      "Intermediate" = "#95a5a6",
      "Strong El Ni\u00f1o" = "#e74c3c"
    )) +
    labs(
      title = "Winter: Temperature vs Rainfall by ENSO Phase",
      subtitle = "Lower-right quadrant = hot and dry (highest agricultural risk)",
      x = "Total Winter Rainfall (mm)", y = "Mean Winter Tmax (\u00b0C)",
      color = "ENSO Phase"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  print(p21)
  ggsave("Q7_enso_compound_winter.png", p21, width = 10, height = 6, dpi = 150)
  
  # Define "bad year" criteria:
  # - Annual rainfall below 25th percentile
  # - Annual mean Tmax above 75th percentile
  # - (For 2000+) Annual mean radiation below 25th percentile
  
  annual_summary <- climate_data %>%
    filter(Ano >= 1917) %>%
    group_by(Ano) %>%
    summarise(
      Total_Rain = sum(Chuva, na.rm = TRUE),
      Mean_TMAX = mean(TMAX, na.rm = TRUE),
      Mean_GlobRad = mean(GlobRad, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate thresholds
  rain_q25 <- quantile(annual_summary$Total_Rain, 0.25, na.rm = TRUE)
  tmax_q75 <- quantile(annual_summary$Mean_TMAX, 0.75, na.rm = TRUE)
  
  annual_summary <- annual_summary %>%
    mutate(
      Low_Rain = Total_Rain < rain_q25,
      High_Temp = Mean_TMAX > tmax_q75,
      Bad_Year = Low_Rain & High_Temp
    )
  
  # Count bad years by decade
  bad_years_decade <- annual_summary %>%
    mutate(Decade = floor(Ano / 10) * 10) %>%
    group_by(Decade) %>%
    summarise(
      N_Years = n(),
      N_Bad = sum(Bad_Year, na.rm = TRUE),
      Pct_Bad = N_Bad / N_Years * 100,
      .groups = "drop"
    )
  
  print(bad_years_decade)
  
  # Plot
  p22 <- ggplot(bad_years_decade, aes(x = factor(Decade), y = Pct_Bad)) +
    geom_col(fill = "#c0392b", alpha = 0.7) +
    labs(
      title = "Probability of a 'Bad Year' by Decade",
      subtitle = "Bad year = low rainfall (< Q25) AND high Tmax (> Q75)",
      x = "Decade", y = "Percentage of Bad Years (%)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p22)
  ggsave("Q7_bad_years.png", p22, width = 10, height = 5, dpi = 150)
  
  # Sugarcane critical period: August-November (tillering to grand growth)
  # Water stress = extended dry spell + high temperature + low humidity
  sugarcane_risk <- climate_data %>%
    filter(Ano >= 1917, Mes %in% c(8, 9, 10, 11)) %>%
    group_by(Ano) %>%
    summarise(
      Mean_TMAX = mean(TMAX, na.rm = TRUE),
      Total_Rain = sum(Chuva, na.rm = TRUE),
      Max_Estiagem = max(Estiagem, na.rm = TRUE),
      Mean_URMED = mean(URMED, na.rm = TRUE),
      # Composite stress index (standardized)
      .groups = "drop"
    ) %>%
    mutate(
      # Z-score standardization for composite index
      Z_Temp = scale(Mean_TMAX)[,1],
      Z_Rain = -scale(Total_Rain)[,1],  # Negative: less rain = more stress
      Z_Drought = scale(Max_Estiagem)[,1],
      Stress_Index = (Z_Temp + Z_Rain + Z_Drought) / 3
    )
  
  # Trend in stress index
  cat("Sugarcane stress index trend:\n")
  print(trend_summary(sugarcane_risk$Ano, sugarcane_risk$Stress_Index))
  
  # Plot
  p23 <- ggplot(sugarcane_risk, aes(x = Ano, y = Stress_Index)) +
    geom_line(color = "grey50", alpha = 0.6) +
    geom_smooth(method = "lm", color = "#c0392b") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "red") +
    annotate("text", x = 1920, y = 1.1, label = "High stress threshold",
             color = "red", size = 3, hjust = 0) +
    labs(
      title = "Sugarcane Water Stress Index (Aug\u2013Nov)",
      subtitle = "Composite of temperature, rainfall deficit, and drought length",
      x = "Year", y = "Stress Index (standardized)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p23)
  ggsave("Q7_sugarcane_stress.png", p23, width = 10, height = 5, dpi = 150)
  
  # ==============================================================================
  # ITEM 8: ANÁLISE DE RISCO DE INCÊNDIO FLORESTAL (METODOLOGIA DE ANGSTRÖM)
  # ==============================================================================
  
  # 1. Cálculo Diário do Índice de Angström e Classificação de Perigo
  # Valores baixos de Angström indicam maior inflamabilidade (Ar seco e quente)
  fire_risk_daily <- climate_data %>%
    filter(Ano >= 1917, !is.na(TMAX), !is.na(URMED)) %>%
    mutate(
      # Cálculo do Índice de Angström
      Angstrom_Index = (URMED / 20) + ((TMAX - 27) / 10),
      
      # Classificação Tradicional do Risco Florestal de Angström
      Classe_Risco = case_when(
        Angstrom_Index > 4.0  ~ "Nulo/Baixo",
        Angstrom_Index > 2.5  ~ "Moderado",
        Angstrom_Index > 2.0  ~ "Alto",
        TRUE                  ~ "Muito Alto a Extremo"
      ),
      Classe_Risco = factor(Classe_Risco, 
                            levels = c("Nulo/Baixo", "Moderado", "Alto", "Muito Alto a Extremo"))
    )
  
  # 2. Agrupamento Anual para Análise de Tendência de Longo Prazo
  fire_risk_annual <- fire_risk_daily %>%
    group_by(Ano) %>%
    summarise(
      Mean_Angstrom = mean(Angstrom_Index, na.rm = TRUE),
      Days_Critical_Risk = sum(Angstrom_Index <= 2.5, na.rm = TRUE), # Dias com risco Alto + Extremo
      Pct_Critical_Days  = (Days_Critical_Risk / n()) * 100,
      .groups = "drop"
    )
  
  # Visualizar as primeiras linhas do output estruturado
  cat("\n--- Primeiras linhas dos indicadores anuais de Risco de Incêndio ---\n")
  print(head(fire_risk_annual))
  
  # 3. Modelagem Linear e Teste de Tendência de Mann-Kendall / Sen's Slope
  lm_angstrom <- lm(Mean_Angstrom ~ Ano, data = fire_risk_annual)
  lm_critical <- lm(Days_Critical_Risk ~ Ano, data = fire_risk_annual)
  
  cat("\n--- Tendência Linear: Índice Médio de Angström (por década) ---\n")
  cat("Variação por década:", coef(lm_angstrom)[2] * 10, "\n")
  
  cat("\n--- Testes de Mann-Kendall não-paramétricos ---\n")
  cat("Índice de Angström Médio:\n")
  print(trend_summary(fire_risk_annual$Ano, fire_risk_annual$Mean_Angstrom))
  
  cat("\nDias Críticos de Risco de Incêndio (Angström <= 2.5):\n")
  print(trend_summary(fire_risk_annual$Ano, fire_risk_annual$Days_Critical_Risk))
  
  # 4. Plotagem Gráfica 8A: Evolução do Número de Dias Críticos por Ano
  p24 <- ggplot(fire_risk_annual, aes(x = Ano, y = Days_Critical_Risk)) +
    geom_line(color = "#d35400", alpha = 0.5) +
    geom_point(color = "#c0392b", alpha = 0.7, size = 1.5) +
    geom_smooth(method = "lm", color = "#2c3e50", fill = "#f8d7da", linewidth = 1.2) +
    labs(
      title = "Evolução do Perigo de Incêndios Florestais em Piracicaba (1917-Presente)",
      subtitle = "Número anual de dias com Risco Alto a Extremo (Índice de Angström \u2264 2.5)",
      x = "Ano", y = "Número de Dias Críticos / Ano"
    ) +
    theme_minimal(base_size = 12)
  
  print(p24)
  ggsave("Q8_dias_criticos_incendio.png", p24, width = 10, height = 5, dpi = 150)
  
  # 5. Plotagem Gráfica 8B: Sazonalidade Histórica e Impacto do ENSO no Risco
  fire_season_enso <- fire_risk_daily %>%
    filter(Ano >= 1950, !is.na(ClassNino)) %>%
    group_by(Ano, Estacao, ClassNino) %>%
    summarise(Dias_Risco = sum(Angstrom_Index <= 2.5, na.rm = TRUE), .groups = "drop")
  
  p25 <- ggplot(fire_season_enso, aes(x = ClassNino, y = Dias_Risco, fill = ClassNino)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.4) +
    facet_wrap(~ Estacao, scales = "free_y") +
    scale_fill_manual(values = c(
      "Strong La Ni\u00f1a" = "#3498db",
      "Intermediate" = "#95a5a6",
      "Strong El Ni\u00f1o" = "#e74c3c"
    )) +
    labs(
      title = "Dias Críticos de Incêndio por Estação e Fase do ENSO",
      subtitle = "Metodologia Angström aplicada ao monitoramento de Clima Florestal",
      x = "Fase do ENSO", y = "Dias com Risco Crítico (\u2264 2.5)", fill = "ENSO"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
  
  print(p25)
  ggsave("Q8_risco_incendio_enso.png", p25, width = 11, height = 7, dpi = 150)
  
  # Teste Estatístico para validação do impacto do ENSO no Risco de Incêndio no Inverno
  cat("\n--- Teste de Kruskal-Wallis: Impacto do ENSO no Risco de Incêndio (Inverno) ---\n")
  winter_fire <- fire_season_enso %>% filter(Estacao == "Winter")
  print(kruskal.test(Dias_Risco ~ ClassNino, data = winter_fire))
  
  