
library(dplyr)
library(stringr)
library(ggplot2)

#### calculcate and plot total effects 

# first: AB key stone predator systems
totaleffectsABpred<-total_effect(fit_semAB_Pred_fulltopdown, n_lags = 4)

# removing unwanted observations 
totaleffectsABpred <- totaleffectsABpred %>%
  filter(lag < 2, 
         total_effect != 0, 
         total_effect != 1.0,
         !(str_detect(to, "_") & str_detect(from, "_") & 
             str_extract(to, "(?<=_).*") != str_extract(from, "(?<=_).*")))

# calculate the mean and standard error (se) by grouping lag, 'to' and 'from' groups
totaleffectsABpred_mean <- totaleffectsABpred %>%
  mutate(
    to_variable = str_extract(to, "^[^_]+"),  
    from_variable = str_sub(from, 1, 8)  
  ) %>%
  group_by(lag, to_variable, from_variable) %>%
summarize(mean_total_effect = mean(total_effect, na.rm = TRUE),
    se_total_effect = sd(total_effect, na.rm = TRUE) / sqrt(n()),
    mean_direct_effect = mean(direct_effect, na.rm = TRUE),
    se_direct_effect = sd(direct_effect, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

# remove observations where salinity is a response 
totaleffectsABpred_mean <- totaleffectsABpred_mean %>%
  filter(to_variable != "Salinity")

# create column for differnet relationship types 
totaleffectsABpred_mean <- totaleffectsABpred_mean %>%
 mutate(relationship_type = case_when(
      from_variable == "PDSI" ~ "PDSI",
      from_variable == "Salinity" ~ "Salinity",
      substr(to_variable, 1, 7) == substr(from_variable, 1, 7) ~ "Density Dependent",
      TRUE ~ "Trophic"))

# calculcate mean absolute values (may help with communication overall variable importance)
totaleffectsABpred_mean <- totaleffectsABpred_mean %>%
  mutate(mean_total_effect_absvalue = abs(mean_total_effect))

relationship_type_averages <- totaleffectsABpred_mean %>%
  group_by(relationship_type) %>%
  summarise(
    avg_abs_effect = mean(mean_total_effect_absvalue, na.rm = TRUE),
    se_abs_effect = sd(mean_total_effect_absvalue, na.rm = TRUE) / sqrt(n())  ) %>%
  arrange(desc(avg_abs_effect))

# creating a color pallette
custom_colors_ABpred <- c(
  "Salinity" = "#f4a261", 
  "PDSI" = "#f28482", 
  "AllMulle" = "#e9c46a", 
  "AllMenha" = "#90a955",
  "BullShar" = "#2a9d8f", 
  "Alligato" = "#264653")

# create vertical breaks for plot
x_breaks <- unique(as.numeric(factor(totaleffectsABpred_mean$to_variable)))
vline_positions <- x_breaks + 0.5  

# create legend labels
lag_labels <- c("0" = "No Lag", "1" = "1 Year Lag")
predictor_labels_ABpred <- c(
  "Alligato" = "Alligator Gar",
  "AllMenha" = "Menhaden",
  "AllMulle" = "Mullet",
  "BullShar" = "Bull Shark")

# turn on to remove prey species
#totaleffectsABpred_mean <- totaleffectsABpred_mean %>%
  #filter(!to_variable %in% c("AllMullet", "AllMenhaden"))

# rename two taxa for consitency
totaleffectsABpred_mean <- totaleffectsABpred_mean %>%
  mutate(to_variable = recode(to_variable,
                              "AllMullet" = "Mullet",
                              "AllMenhaden" = "Menhaden"))

ABpred_totaleffectsplot<-ggplot(totaleffectsABpred_mean, aes(x = to_variable, y = mean_total_effect, fill = from_variable)) +
  geom_errorbar(aes(
    ymin = mean_total_effect - se_total_effect,
    ymax = mean_total_effect + se_total_effect 
  ), width = 0.25, position = position_dodge(0.9)) +  
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7, color = "black") +
  geom_vline(xintercept = vline_positions, linetype = "dashed", color = "gray50") +
  facet_wrap(~ lag, labeller = labeller(lag = lag_labels)) +
  labs(title = "Mission-Aransas Estuary - Keystone Predator System",
       y = "Mean Total Effect", x = "Response", fill = "Predictor") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, vjust = 0.6),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)) +
  scale_fill_manual(values = custom_colors_ABpred, labels = predictor_labels_ABpred)

ggsave("ABpred_totaleffectsplot.png", ABpred_totaleffectsplot, dpi = 150, bg = "white",
       width = 2000, height = 1000, units = "px") 

########################



# second: AB sciaenid  system
totaleffectsABSciaenid<-total_effect(fit_semAB_Sciaenid_notrophics, n_lags = 4)

# removing unwanted observations 
totaleffectsABSciaenid <- totaleffectsABSciaenid %>%
  filter(lag < 2, 
         total_effect != 0, 
         total_effect != 1.0,
         !(str_detect(to, "_") & str_detect(from, "_") & 
             str_extract(to, "(?<=_).*") != str_extract(from, "(?<=_).*")))

# calculate the mean and standard error (se) by grouping lag, 'to' and 'from' groups
totaleffectsABSciaenid_mean <- totaleffectsABSciaenid %>%
  mutate(
    to_variable = str_extract(to, "^[^_]+"),  
    from_variable = str_sub(from, 1, 8)  
  ) %>%
  group_by(lag, to_variable, from_variable) %>%
  summarize(mean_total_effect = mean(total_effect, na.rm = TRUE),
            se_total_effect = sd(total_effect, na.rm = TRUE) / sqrt(n()),
            mean_direct_effect = mean(direct_effect, na.rm = TRUE),
            se_direct_effect = sd(direct_effect, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# remove observations where salinity is a response 
totaleffectsABSciaenid_mean <- totaleffectsABSciaenid_mean %>%
  filter(to_variable != "Salinity")

# create column for different relationship types 
totaleffectsABSciaenid_mean <- totaleffectsABSciaenid_mean %>%
  mutate(relationship_type = case_when(
    from_variable == "PDSI" ~ "PDSI",
    from_variable == "Salinity" ~ "Salinity",
    substr(to_variable, 1, 7) == substr(from_variable, 1, 7) ~ "Density Dependent",
    TRUE ~ "Trophic"))

# calculcate mean absolute values (may help with communication overall variable importance)
totaleffectsABSciaenid_mean <- totaleffectsABSciaenid_mean %>%
  mutate(mean_total_effect_absvalue = abs(mean_total_effect))

relationship_type_averages <- totaleffectsABSciaenid_mean %>%
  group_by(relationship_type) %>%
  summarise(
    avg_abs_effect = mean(mean_total_effect_absvalue, na.rm = TRUE),
    se_abs_effect = sd(mean_total_effect_absvalue, na.rm = TRUE) / sqrt(n())  ) %>%
  arrange(desc(avg_abs_effect))

# creating a color pallette
custom_colors_ABSciaenid <- c(
  "Salinity" = "#f4a261", 
  "PDSI" = "#f28482", 
  "RedDrum_" = "#2a9d8f", 
  "SpottedS" = "#264653")

# create vertical breaks for plot
x_breaks <- unique(as.numeric(factor(totaleffectsABSciaenid_mean$to_variable)))
vline_positions <- x_breaks + 0.5  

# create legend labels
lag_labels <- c("0" = "No Lag", "1" = "1 Year Lag")
predictor_labels_ABSciaenid <- c(
  "RedDrum_" = "Red Drum",
  "SpottedS" = "Spotted Seatrout")

# turn on to remove prey species
#totaleffectsABSciaenid_mean<- totaleffectsABSciaenid_mean%>%
  #filter(!to_variable %in% c("Atlanticcroaker", "BlueCrabSmall"))

ABSciaenid_totaleffectsplot<-ggplot(totaleffectsABSciaenid_mean, aes(x = to_variable, y = mean_total_effect, fill = from_variable)) +
  geom_errorbar(aes(
    ymin = mean_total_effect - se_total_effect,
    ymax = mean_total_effect + se_total_effect 
  ), width = 0.25, position = position_dodge(0.9)) +  
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7, color = "black") +
  geom_vline(xintercept = vline_positions, linetype = "dashed", color = "gray50") +
  facet_wrap(~ lag, labeller = labeller(lag = lag_labels)) +
  labs(title = "Fina - Sciaenid System",
       y = "Mean Total Effect", x = "Response", fill = "Predictor") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, vjust = 0.6),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = custom_colors_ABSciaenid, labels = predictor_labels_ABSciaenid)+
  scale_x_discrete(labels = c(
    "Atlanticcroaker" = "AtlanticCroaker",
    "BlueCrabSmall" = "BlueCrab"))

ggsave("ABSciaenid_totaleffectsplot.png", ABSciaenid_totaleffectsplot, dpi = 150, bg = "white",
       width = 2000, height = 1000, units = "px") 

##########################


# second: GB sciaenid system
totaleffectsGBSciaenid<-total_effect(fit_semGB_Sciaenid_fulltopdown, n_lags = 4)

# removing unwanted observations 
totaleffectsGBSciaenid <- totaleffectsGBSciaenid %>%
  filter(lag < 2, 
         total_effect != 0, 
         total_effect != 1.0,
         !(str_detect(to, "_") & str_detect(from, "_") & 
             str_extract(to, "(?<=_).*") != str_extract(from, "(?<=_).*")))

# calculate the mean and standard error (se) by grouping lag, 'to' and 'from' groups
totaleffectsGBSciaenid_mean <- totaleffectsGBSciaenid %>%
  mutate(
    to_variable = str_extract(to, "^[^_]+"),  
    from_variable = str_sub(from, 1, 8)  
  ) %>%
  group_by(lag, to_variable, from_variable) %>%
  summarize(mean_total_effect = mean(total_effect, na.rm = TRUE),
            se_total_effect = sd(total_effect, na.rm = TRUE) / sqrt(n()),
            mean_direct_effect = mean(direct_effect, na.rm = TRUE),
            se_direct_effect = sd(direct_effect, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# remove observations where salinity is a response 
totaleffectsGBSciaenid_mean <- totaleffectsGBSciaenid_mean %>%
  filter(to_variable != "Salinity")

# create column for different relationship types 
totaleffectsGBSciaenid_mean <- totaleffectsGBSciaenid_mean %>%
  mutate(relationship_type = case_when(
    from_variable == "PDSI" ~ "PDSI",
    from_variable == "Salinity" ~ "Salinity",
    substr(to_variable, 1, 7) == substr(from_variable, 1, 7) ~ "Density Dependent",
    TRUE ~ "Trophic"))

# calculcate mean absolute values (may help with communication overall variable importance)
totaleffectsGBSciaenid_mean <- totaleffectsGBSciaenid_mean %>%
  mutate(mean_total_effect_absvalue = abs(mean_total_effect))

relationship_type_averages <- totaleffectsGBSciaenid_mean %>%
  group_by(relationship_type) %>%
  summarise(
    avg_abs_effect = mean(mean_total_effect_absvalue, na.rm = TRUE),
    se_abs_effect = sd(mean_total_effect_absvalue, na.rm = TRUE) / sqrt(n())  ) %>%
  arrange(desc(avg_abs_effect))

# creating a color pallette
custom_colors_GBSciaenid <- c(
  "Salinity" = "#f4a261", 
  "PDSI" = "#f28482", 
  "Atlantic" = "#e9c46a", 
  "BlueCrab" = "#90a955",
  "RedDrum_" = "#2a9d8f", 
  "SpottedS" = "#264653")

# create vertical breaks for plot
x_breaks <- unique(as.numeric(factor(totaleffectsGBSciaenid_mean$to_variable)))
vline_positions <- x_breaks + 0.5  

# create legend labels
lag_labels <- c("0" = "No Lag", "1" = "1 Year Lag")
predictor_labels_GBSciaenid <- c(
  "RedDrum_" = "Red Drum",
  "Atlantic" = "Atlantic Croaker",
  "BlueCrab" = "Blue Crab",
  "SpottedS" = "Spotted Seatrout")

# turn on to remove prey species
#totaleffectsGBSciaenid_mean<- totaleffectsGBSciaenid_mean%>%
  #filter(!to_variable %in% c("Atlanticcroaker", "BlueCrabSmall"))

GBSciaenid_totaleffectsplot<-ggplot(totaleffectsGBSciaenid_mean, aes(x = to_variable, y = mean_total_effect, fill = from_variable)) +
  geom_errorbar(aes(
    ymin = mean_total_effect - se_total_effect,
    ymax = mean_total_effect + se_total_effect 
  ), width = 0.25, position = position_dodge(0.9)) +  
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7, color = "black") +
  geom_vline(xintercept = vline_positions, linetype = "dashed", color = "gray50") +
  facet_wrap(~ lag, labeller = labeller(lag = lag_labels)) +
  labs(title = "Trinity-San Jacinto Estuary - Sciaenid System",
       y = "Mean Total Effect", x = "Response", fill = "Predictor") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, vjust = 0.6),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = custom_colors_GBSciaenid, labels = predictor_labels_GBSciaenid)+
  scale_x_discrete(labels = c(
    "Atlanticcroaker" = "AtlanticCroaker",
    "BlueCrabSmall" = "BlueCrab"))

ggsave("GBSciaenid_totaleffectsplot.png", GBSciaenid_totaleffectsplot, dpi = 150, bg = "white",
       width = 2000, height = 1000, units = "px") 

##########################



# last: GB key stone predator systems
totaleffectsGBpred<-total_effect(fit_semGB_Pred_fullbottomup, n_lags = 4)

# removing unwanted observations 
totaleffectsGBpred <- totaleffectsGBpred %>%
  filter(lag < 2, 
         total_effect != 0, 
         total_effect != 1.0,
         !(str_detect(to, "_") & str_detect(from, "_") & 
             str_extract(to, "(?<=_).*") != str_extract(from, "(?<=_).*")))

# calculate the mean and standard error (se) by grouping lag, 'to' and 'from' groups
totaleffectsGBpred_mean <- totaleffectsGBpred %>%
  mutate(
    to_variable = str_extract(to, "^[^_]+"),  
    from_variable = str_sub(from, 1, 3)  
  ) %>%
  group_by(lag, to_variable, from_variable) %>%
  summarize(mean_total_effect = mean(total_effect, na.rm = TRUE),
            se_total_effect = sd(total_effect, na.rm = TRUE) / sqrt(n()),
            mean_direct_effect = mean(direct_effect, na.rm = TRUE),
            se_direct_effect = sd(direct_effect, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# remove observations where salinity is a response 
totaleffectsGBpred_mean <- totaleffectsGBpred_mean %>%
  filter(to_variable != "Salinity")

# create column for different relationship types 
totaleffectsGBpred_mean <- totaleffectsGBpred_mean %>%
  mutate(relationship_type = case_when(
    from_variable == "PDS" ~ "PDSI",
    from_variable == "Sal" ~ "Sal",
    substr(to_variable, 1, 3) == substr(from_variable, 1, 3) ~ "Density Dependent",
    TRUE ~ "Trophic"))

# calculcate mean absolute values (may help with communication overall variable importance)
totaleffectsGBpred_mean <- totaleffectsGBpred_mean %>%
  mutate(mean_total_effect_absvalue = abs(mean_total_effect))

relationship_type_averages <- totaleffectsGBpred_mean %>%
  group_by(relationship_type) %>%
  summarise(
    avg_abs_effect = mean(mean_total_effect_absvalue, na.rm = TRUE),
    se_abs_effect = sd(mean_total_effect_absvalue, na.rm = TRUE) / sqrt(n())  ) %>%
  arrange(desc(avg_abs_effect))

# creating a color pallette
custom_colors_GBpred <- c(
  "Sal" = "#f4a261", 
  "PDS" = "#f28482", 
  "Mul" = "#e9c46a", 
  "Men" = "#90a955",
  "Bul" = "#2a9d8f", 
  "All" = "#264653")

# create vertical breaks for plot
x_breaks <- unique(as.numeric(factor(totaleffectsGBpred_mean$to_variable)))
vline_positions <- x_breaks + 0.5  

# create legend labels
lag_labels <- c("0" = "No Lag", "1" = "1 Year Lag")
predictor_labels_GBpred <- c(
  "Sal" = "Salinity",
  "PDS" = "PDSI",
  "All" = "Alligator Gar",
  "Men" = "Menhaden",
  "Mul" = "Mullet",
  "Bul" = "Bull Shark")

# turn on to remove prey species
#totaleffectsGBpred_mean <- totaleffectsGBpred_mean %>%
  #filter(!to_variable %in% c("Mullet", "Menhaden"))

GBpred_totaleffectsplot<-ggplot(totaleffectsGBpred_mean, aes(x = to_variable, y = mean_total_effect, fill = from_variable)) +
  geom_errorbar(aes(
    ymin = mean_total_effect - se_total_effect,
    ymax = mean_total_effect + se_total_effect 
  ), width = 0.25, position = position_dodge(0.9)) +  
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7, color = "black") +
  geom_vline(xintercept = vline_positions, linetype = "dashed", color = "gray50") +
  facet_wrap(~ lag, labeller = labeller(lag = lag_labels)) +
  labs(title = "Trinity-San Jacinto Estuary - Keystone Predator System",
  y = "Mean Total Effect", x = "Response", fill = "Predictor") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, vjust = 0.6),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = custom_colors_GBpred, labels = predictor_labels_GBpred)

ggsave("GBpred_totaleffectsplot.png", GBpred_totaleffectsplot, dpi = 150, bg = "white",
       width = 2000, height = 1000, units = "px") 

