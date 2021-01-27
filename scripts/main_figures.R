library(tidyverse)
library(ggeffects)
library(ggpubr)
library(here)
library(sjPlot)
library(lme4)
library(SplinesUtils)
library(cowplot)


#load data
rawdata <- readr::read_csv(here::here("data", "raw_data.csv"))
cleandata <- readr::read_csv(here::here("data", "cleaned_data.csv"))
manualdata <- readr::read_csv(here::here("data", "manual_data.csv"))


############################################
################  Figure 2  ################
############################################          

# Quantitative measurements of animal size (raw data)

raw.length <- rawdata %>%
  dplyr::filter(!replicate %in% c("no food", "no worms")) %>%
  dplyr::mutate(hour = as.numeric(timepoint)) %>%
  ggplot() +
  aes(x = hour, y = TOF) +
  geom_jitter(size = 0.01, width = 0.2, alpha = 0.01)+
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  labs(x="Time (hours)", y = expression(paste("Length (", mu,"m)"))) +
  theme_cowplot(font_size = 12, rel_small = 10/12)

raw.width <- rawdata %>%
  dplyr::filter(!replicate %in% c("no food", "no worms")) %>%
  dplyr::mutate(hour = as.numeric(timepoint)) %>%
  ggplot() +
  aes(x = hour, y = norm.EXT) +
  geom_jitter(size = 0.01, width = 0.2, alpha = 0.01)+
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  labs(x="Time (hours)", y = expression(paste("Width (", mu,"m)"))) +
  theme_cowplot(font_size = 12, rel_small = 10/12)

raw.volume <- rawdata %>%
  dplyr::filter(!replicate %in% c("no food", "no worms")) %>%
  dplyr::mutate(hour = as.numeric(timepoint),
                volume = pi*(norm.EXT/2)^2*TOF) %>%
  ggplot() +
  aes(x = hour, y = volume) +
  geom_jitter(size = 0.01, width = 0.2, alpha = 0.01)+
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  labs(x="Time (hours)", y = expression(paste("log Volume (", mu,m^3,")"))) +
  theme_cowplot(font_size = 12, rel_small = 10/12) +
  scale_y_log10()

clean.length <- cleandata %>%
  dplyr::mutate(hour = as.numeric(timepoint)) %>%
  ggplot() +
  aes(x = hour, y = TOF) +
  geom_jitter(size = 0.01, width = 0.2, alpha = 0.01)+
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  labs(x="Time (hours)", y = expression(paste("Length (", mu,"m)"))) +
  theme_cowplot(font_size = 12, rel_small = 10/12)

clean.width <- cleandata %>%
  dplyr::mutate(hour = as.numeric(timepoint)) %>%
  ggplot() +
  aes(x = hour, y = norm.EXT) +
  geom_jitter(size = 0.01, width = 0.2, alpha = 0.01)+
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  labs(x="Time (hours)", y = expression(paste("Width (", mu,"m)"))) +
  theme_cowplot(font_size = 12, rel_small = 10/12) 

clean.volume <- cleandata %>%
  dplyr::mutate(hour = as.numeric(timepoint),
                volume = pi*(norm.EXT/2)^2*TOF) %>%
  ggplot() +
  aes(x = hour, y = volume) +
  geom_jitter(size = 0.01, width = 0.2, alpha = 0.01)+
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  labs(x="Time (hours)", y = expression(paste("log Volume (", mu,m^3,")"))) +
  theme_cowplot(font_size = 12, rel_small = 10/12) +
  scale_y_log10()

fig2 <- cowplot::plot_grid(raw.length, clean.length, raw.width, clean.width, raw.volume, clean.volume, 
                           nrow = 3, ncol = 2, align = "hv", labels = c("A", "D", "B", "E", "C", "F"), label_size = 12)
ggsave(here::here("main_figures", "figure2.jpeg"), plot = fig2, device = "jpeg", width = 25, height = 20, units = "cm", dpi = 300)

############################################
################  Figure 3  ################
############################################ 

# local mins + inflection pts

f.allmin_max <- function() {
  #step 1
  ## calculating residuals
  res <- cleandata %>%
    dplyr::mutate(hour = as.numeric(timepoint), norm.red = red/EXT, well = as.factor(col)) %>%
    aov(norm.red ~ replicate, data = .) %>% residuals()
  ## summarizing data 
  medianplate <- cleandata %>%
    dplyr::mutate(hour = as.numeric(timepoint), volume = pi*(norm.EXT/2)^2*TOF, norm.red = red/EXT,
                  resids = res) %>%
    dplyr::group_by(replicate, hour, timepoint) %>%
    dplyr::summarize(median.residual.red = median(resids),
                     median.norm.red = median(norm.red),
                     median.Length = median(TOF), 
                     median.Width = median(norm.EXT),
                     median.Volume = median(volume), .groups = "drop")

  #step 2: kernel regression smoothing using a local plug-in bandwidth
  step1_deriv <- lokern::lokerns(x=medianplate$hour, y=medianplate$median.residual.red, deriv=1, hetero = T)
  #step 3: backsolve to find where derivative = 0
  step2_backsolve <- CubicInterpSplineAsPiecePoly(x=step1_deriv$x.out,y=step1_deriv$est,method="natural")
  time_min.max <- dplyr::tibble(Hour = solve(step2_backsolve,b=0)) %>%
    dplyr::filter(Hour > 12, Hour < 55)
  #step 4: 
  step0 <- lokern::lokerns(x=medianplate$hour, y=medianplate$median.residual.red, deriv=0, hetero = T)
  plot(step0); abline(v=time_min.max$Hour, col="blue")
  
   #put it all together
  values <- time_min.max %>%
    dplyr::mutate(Fluorescence = c("min","max","min","max","min","max","min","max"),
                  Peak = c(1,2,3,4,5,6,7,8)) %>%
    dplyr::select(Peak, Fluorescence, Hour)
  
  return(values)
}

min_max <- f.allmin_max()

f.allinflection_pts <- function() {
  #step 1
  ## calculating residuals
  res <- cleandata %>%
    dplyr::mutate(hour = as.numeric(timepoint), norm.red = red/EXT, well = as.factor(col)) %>%
    aov(norm.red ~ replicate, data = .) %>% residuals()
  ## summarizing data 
  medianplate <- cleandata %>%
    dplyr::mutate(hour = as.numeric(timepoint), volume = pi*(norm.EXT/2)^2*TOF, norm.red = red/EXT,
                  resids = res) %>%
    dplyr::group_by(replicate, hour, timepoint) %>%
    dplyr::summarize(median.residual.red = median(resids),
                     median.norm.red = median(norm.red),
                     median.Length = median(TOF), 
                     median.Width = median(norm.EXT),
                     median.Volume = median(volume), .groups = "drop")
  
  #step 2: kernel regression smoothing using a local plug-in bandwidth
  step1_deriv <- lokern::lokerns(x=medianplate$hour, y=medianplate$median.residual.red, deriv=2, hetero = T)
  #step 3: backsolve to find where derivative = 0
  step2_backsolve <- CubicInterpSplineAsPiecePoly(x=step1_deriv$x.out,y=step1_deriv$est,method="natural")
  time_min.max <- dplyr::tibble(Hour = solve(step2_backsolve,b=0)) %>%
    dplyr::filter(Hour > 13, Hour < 58)
  #step 4: 
  step0 <- lokern::lokerns(x=medianplate$hour, y=medianplate$median.residual.red, deriv=0, hetero = T)
  plot(step0); abline(v=time_min.max$Hour, col="blue")
  
  #put it all together
  values <- time_min.max %>%
    dplyr::mutate(Fluorescence = c("min","max","min","max","min","max","min","max"),
                  Peak = c(1,2,3,4,5,6,7,8)) %>%
    dplyr::select(Peak, Fluorescence, Hour)
  
  return(values)
}

inflection <- f.allinflection_pts()

# plot
mins <- min_max %>% dplyr::filter(Fluorescence == "min") %>% pull(Hour)
inflection_min <- inflection %>% dplyr::filter(Fluorescence == "min") %>% pull(Hour)
inflection_max <- inflection %>% dplyr::filter(Fluorescence == "max") %>% pull(Hour)
ymin <- c(-Inf, -Inf, -Inf, -Inf); ymax <- c(Inf, Inf, Inf, Inf)
plt <- data.frame(inflection_min, inflection_max, ymin, ymax)

red <- cleandata %>%
  dplyr::group_by(replicate, timepoint, row, col) %>%
  dplyr::summarize(median.red = median(red/EXT), .groups = "drop") %>%
  ggplot(.) +
  aes(x = as.numeric(timepoint), y = median.red) + 
  geom_boxplot(aes(group = timepoint), outlier.alpha = 0) +
  geom_jitter(size=0.2, width=0.2, alpha=0.3) +
  geom_smooth(span = 0.1, se = F, size = 0.8, method = "loess") +
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  geom_vline(xintercept=c(mins), color="red", size = 0.3) +
  #geom_vline(xintercept=c(maxs), color="blue", size = 0.3, alpha = 0.5) +
  #geom_rect(data = plt, aes(xmin=inflection_min, xmax=inflection_max, ymin = ymin, ymax = ymax), fill = "gray", alpha = 0.2, inherit.aes = F) +
  theme_cowplot(font_size = 12, rel_small = 10/12) +
  labs(x="Time (h)", y = "Median Red")

length <- cleandata %>%
  dplyr::group_by(timepoint, replicate, row, col) %>%
  dplyr::summarize(median.Length = median(TOF), .groups = "drop") %>%
  ggplot(., aes(x = as.numeric(timepoint), y = median.Length)) +
  geom_boxplot(aes(group = timepoint), outlier.alpha = 0) +
  geom_jitter(size=0.1, width=0.2, alpha=0.3) +
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  theme_cowplot(font_size = 12, rel_small = 10/12) +
  geom_vline(xintercept=c(mins), color="red", size = 0.3) +
  #geom_rect(data = plt, aes(xmin=inflection_min, xmax=inflection_max, ymin = ymin, ymax = ymax), fill = "gray", alpha = 0.2, inherit.aes = F) +
  labs(x="Time (h)", y = expression(paste("Median Length (", mu,"m)")))

width <- cleandata %>%
  dplyr::group_by(timepoint, replicate, row, col) %>%
  dplyr::summarize(median.Width = median(norm.EXT), .groups = "drop") %>%
  ggplot(.) +
  aes(x = as.numeric(timepoint), y = median.Width) +
  geom_boxplot(aes(group = timepoint), outlier.alpha = 0) +
  geom_jitter(size=0.1, width=0.2, alpha=0.3) +
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  theme_cowplot(font_size = 12, rel_small = 10/12) +
  geom_vline(xintercept=c(mins), color="red", size = 0.3) +
  #geom_rect(data = plt, aes(xmin=inflection_min, xmax=inflection_max, ymin = ymin, ymax = ymax), fill = "gray", alpha = 0.2, inherit.aes = F) +
  labs(x="Time (h)", y = expression(paste("Median Width (", mu,"m)")))

volume <- cleandata %>%
  dplyr::mutate(volume = pi*(norm.EXT/2)^2*TOF) %>%
  dplyr::group_by(timepoint, replicate, row, col) %>%
  dplyr::summarize(median.Volume = median(volume),
                   median.logVolume = median(log(volume)), .groups = "drop") %>%
  ggplot(.) +
  aes(x = as.numeric(timepoint), y = median.logVolume) +
  geom_boxplot(aes(group = timepoint), outlier.alpha = 0) +
  geom_jitter(size=0.1, width=0.2, alpha=0.3) +
  scale_x_continuous(breaks = seq(0, 70, 5)) +
  theme_cowplot(font_size = 12, rel_small = 10/12) +
  geom_vline(xintercept=c(mins), color="red", size = 0.3) +
  #scale_y_log10() +
  #geom_rect(data = plt, aes(xmin=inflection_min, xmax=inflection_max, ymin = ymin, ymax = ymax), fill = "gray", alpha = 0.2, inherit.aes = F) +
  labs(x="Time (h)", y=expression(paste("Median log Volume (", mu, m^3,")")))

fig3 <- cowplot::plot_grid(red, length, width, volume, labels = c("A", "B", "C", "D"), ncol = 1, nrow = 4, align = "hv", label_size = 12)

ggsave(here::here("main_figures", "figure3.jpeg"), plot = fig3, device = "jpeg", width = 20, height = 25, units = "cm", dpi = 300)


############################################
################  Figure 5  ################
############################################ 

### molt density plot
library(ggridges)
cleandata %>%
  dplyr::mutate(hour = as.numeric(timepoint),
                norm.red.area = red/EXT, volume = pi*(norm.EXT/2)^2*TOF,
                Quiescent = dplyr::case_when(hour > 13 & norm.red.area <= 0.05 ~ "yes",
                                             TRUE ~ "no")) %>%
  dplyr::filter(timepoint >= 32, timepoint <= 41) %>%
  ggplot() +
  aes(x = TOF, y = timepoint, group = timepoint, point_color = Quiescent) +
  geom_density_ridges(alpha = 0.2, quantile_lines = T, quantiles = 5, jittered_points = T, 
                      position = position_raincloud(seed = 2, height = 0.3), 
                      point_alpha = 0.4, point_size = 0.1, scale = 0.55) +
  scale_discrete_manual("point_color", values = c(yes = "red", no = "black")) +
  #scale_point_color_gradient(low = "yellow", high = "darkblue", guide = "legend") +
  theme_cowplot() + 
  #geom_hline(yintercept = c("36", "48"), color = "red", size = 0.5) +
  labs(y ="Time (h)", x = expression(paste("Length")))



cleandata %>%
  dplyr::filter(timepoint >= 32, timepoint <= 40) %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(volume = pi*(norm.EXT/2)^2*TOF,
                q20_TOF = as.numeric(stats::quantile(TOF, probs = 0.2, na.rm = TRUE)[1]),
                q40_TOF = as.numeric(stats::quantile(TOF, probs = 0.4, na.rm = TRUE)[1]),
                q60_TOF = as.numeric(stats::quantile(TOF, probs = 0.6, na.rm = TRUE)[1]),
                q80_TOF = as.numeric(stats::quantile(TOF, probs = 0.8, na.rm = TRUE)[1])) %>%
  dplyr::mutate(quant = dplyr::case_when(TOF <= q20_TOF ~ "Q1",
                                         TOF > q20_TOF & TOF <= q40_TOF ~ "Q2",
                                         TOF > q40_TOF & TOF <= q60_TOF ~ "Q3",
                                         TOF > q60_TOF & TOF <= q80_TOF ~ "Q4",
                                         TOF >= q80_TOF ~ "Q5")) %>%
  dplyr::group_by(timepoint, quant) %>%
  dplyr::mutate(hour = as.numeric(timepoint),
                norm.red.area = red/EXT,
                Quiescent = dplyr::case_when(norm.red.area <= 0.05 ~ "yes",
                                             TRUE ~ "no"),
                feeding = sum(Quiescent == "no"),
                notfeeding = sum(Quiescent == "yes"),
                `% Quiescent` = (notfeeding/(notfeeding + feeding))*100,
                count = n()) %>%
  dplyr::ungroup() %>% 
  ggplot() +
  aes(x = TOF, y = factor(hour), group = timepoint, point_color = `% Quiescent`) +
  geom_density_ridges(alpha = 0.2, quantile_lines = T, quantiles = 5, jittered_points = T, 
                      position = position_points_sina(seed = 2),
                      point_alpha = 0.4, point_size = 0.4, scale = 0.99) +
  scale_point_color_gradient(name = "% Quiescent", low = "black", high = "red",
                             breaks = c(10,20,40,60,80), 
                             labels = c("10","20","40","60","80")) +
  theme_cowplot() + 
  guides(point_color = guide_legend(override.aes = list(size = 0.001, point_size = 2, point_alpha = 1))) +
  labs(y ="Time (h)", x = expression(paste("Length")))








