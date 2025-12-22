
# custom function to get predictions during the 3 driest and the 3 wettest years (according to PDSI)
process_PDSI_data <- function(df, export_name) {
  library(dplyr)
  
  # Step 1: Filter the dataframe to keep only rows where Var2 is 'PDSI'
  df_PDSI <- df %>%
    filter(Var2 == "PDSI") %>%
    arrange(desc(obs))  # Sort by observed PDSI values in descending order
  
  # Step 2: Identify the 3 years with the highest and lowest observed PDSI values
  highest_years <- head(df_PDSI, 3)  # Top 3 highest PDSI years
  lowest_years <- tail(df_PDSI, 3)   # Bottom 3 lowest PDSI years
  
  # Get the years for the highest and lowest PDSI
  highest_years_list <- unique(highest_years$Var1)
  lowest_years_list <- unique(lowest_years$Var1)
  
  # Combine the highest and lowest years
  years_to_keep <- c(highest_years_list, lowest_years_list)
  
  # Step 3: Filter the original dataframe to keep rows for these years across all variables
  df_PDSI_filtered <- df %>%
    filter(Var1 %in% years_to_keep)
  
  top_bottom_years <- df_PDSI_filtered %>%
    filter(Var2 == "PDSI") %>%
    select(Var1, obs) %>%
    arrange(desc(obs)) %>%
    slice(c(1:3, (n() - 2):n()))  # Get the 3 highest and 3 lowest PDSI years
  
  # Step 4: Create the 'weather' column based on the PDSI values of those years
  df_PDSI_filtered <- df_PDSI_filtered %>%
    mutate(weather = ifelse(Var1 %in% top_bottom_years$Var1[top_bottom_years$obs < 0], "Severely Dry", 
                            ifelse(Var1 %in% top_bottom_years$Var1[top_bottom_years$obs > 0], "Severely Wet", NA)))
  
  # Step 5: Calculate mean and standard error of 'est' and 'obs' for each species/minor bay and weather condition
  df_summary <- df_PDSI_filtered %>%
    filter(Var2 != "PDSI") %>%  # Exclude the PDSI rows from the summary
    group_by(Var2, weather) %>%
    summarise(
      mean_est = mean(est, na.rm = TRUE),
      se_est = sd(est, na.rm = TRUE) / sqrt(n()),
      mean_obs = mean(obs, na.rm = TRUE),
      se_obs = sd(obs, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup()
  
  # Step 6: Identify all species and minor bays dynamically
  excluded_vars <- c("PDSI")  # Exclude PDSI from species list
  all_species <- unique(df$Var2[!df$Var2 %in% excluded_vars])  # Identify species
  minor_bays <- unique(sub(".*_", "", all_species))  # Extract minor bay names
  
  # Step 7: Get mean PDSI values for Severely Dry and Severely Wet years
  mean_PDSI_dry <- mean(df_PDSI_filtered$obs[df_PDSI_filtered$Var2 == "PDSI" & df_PDSI_filtered$weather == "Severely Dry"], na.rm = TRUE)
  mean_PDSI_wet <- mean(df_PDSI_filtered$obs[df_PDSI_filtered$Var2 == "PDSI" & df_PDSI_filtered$weather == "Severely Wet"], na.rm = TRUE)
  
  # Step 8: Assign the summary tibble to the user-specified name
  assign(export_name, df_summary, envir = .GlobalEnv)
  
  # Return a list of results
  return(list(
    summary = df_summary,
    mean_PDSI_dry = mean_PDSI_dry,
    mean_PDSI_wet = mean_PDSI_wet,
    species = all_species,
    minor_bays = minor_bays
  ))
}

# Store the result of the function
results_AB_pred <- process_PDSI_data(loodf_AB_Pred, "loodf_AB_Pred_PDSI_forplotting")
results_GB_pred <- process_PDSI_data(loodf_GB_Pred, "loodf_GB_Pred_PDSI_forplotting")
results_AB_Sciaenid <- process_PDSI_data(loodf_AB_Sciaenid, "loodf_AB_Sciaenid_PDSI_forplotting")
results_GB_Sciaenid <- process_PDSI_data(loodf_GB_Sciaenid, "loodf_GB_Sciaenid_PDSI_forplotting")


# Get mean PDSI values for Severely Dry and Severely Wet years
mean_PDSI_dry_AB <- results_AB_pred$mean_PDSI_dry
mean_PDSI_wet_AB <- results_AB_pred$mean_PDSI_wet
mean_PDSI_dry_GB <- results_GB_pred$mean_PDSI_dry
mean_PDSI_wet_GB <- results_GB_pred$mean_PDSI_wet

# Remove and rename rows
loodf_AB_Pred_PDSI_forplotting_filtered  <- loodf_AB_Pred_PDSI_forplotting%>%
  filter(!Var2 %in% c("PDSI", "Salinity_CopanoBay", "Salinity_MesquiteBay", "Salinity_AransasBay")) %>%
  mutate(
    Var2 = recode(Var2,
                  "BullShark_AransasBay" = "Bull Shark (AB)",
                  "AllMenhaden_AransasBay" = "Menhaden (AB)",
                  "AlligatorGar_AransasBay" = "Alligator Gar (AB)",
                  "AllMullet_AransasBay" = "Mullet (AB)",
                  "BullShark_CopanoBay" = "Bull Shark (CB)",
                  "AllMenhaden_CopanoBay" = "Menhaden (CB)",
                  "AlligatorGar_CopanoBay" = "Alligator Gar (CB)",
                  "AllMullet_CopanoBay" = "Mullet (CB)",
                  "BullShark_MesquiteBay" = "Bull Shark (MB)",
                  "AllMenhaden_MesquiteBay" = "Menhaden (MB)",
                  "AlligatorGar_MesquiteBay" = "Alligator Gar (MB)",
                  "AllMullet_MesquiteBay" = "Mullet (MB)"
    )
  )


loodf_GB_Pred_PDSI_forplotting_filtered  <- loodf_GB_Pred_PDSI_forplotting%>%
  filter(!Var2 %in% c("PDSI", "Salinity_WestBay", "Salinity_EastBay", "Salinity_TrinityBay", "Salinity_GalvestonBay")) %>%
  mutate(
    Var2 = recode(Var2,
                  "BullShark_WestBay" = "Bull Shark (WB)",
                  "Menhaden_WestBay" = "Menhaden (WB)",
                  "AlligatorGar_WestBay" = "Alligator Gar (WB)",
                  "Mullet_WestBay" = "Mullet (WB)",
                  "BullShark_EastBay" = "Bull Shark (EB)",
                  "Menhaden_EastBay" = "Menhaden (EB)",
                  "AlligatorGar_EastBay" = "Alligator Gar (EB)",
                  "Mullet_EastBay" = "Mullet (EB)",
                  "BullShark_TrinityBay" = "Bull Shark (TB)",
                  "Menhaden_TrinityBay" = "Menhaden (TB)",
                  "AlligatorGar_TrinityBay" = "Alligator Gar (TB)",
                  "Mullet_TrinityBay" = "Mullet (TB)",
                  "BullShark_GalvestonBay" = "Bull Shark (GB)",
                  "Menhaden_GalvestonBay" = "Menhaden (GB)",
                  "AlligatorGar_GalvestonBay" = "Alligator Gar (GB)",
                  "Mullet_GalvestonBay" = "Mullet (GB)"
    )
  )

loodf_GB_Sciaenid_PDSI_forplotting_filtered  <- loodf_GB_Sciaenid_PDSI_forplotting%>%
  filter(!Var2 %in% c("PDSI", "Salinity_WestBay", "Salinity_EastBay", "Salinity_TrinityBay", "Salinity_GalvestonBay")) %>%
  mutate(
    Var2 = recode(Var2,
                  "RedDrum_WestBay" = "Red Drum (WB)",
                  "BlueCrabSmall_WestBay" = "Blue Crab  (WB)",
                  "SpottedSeatrout_WestBay" = "Spotted Seatrout (WB)",
                  "Atlanticcroaker_WestBay" = "Atlantic Croaker (WB)",
                  "RedDrum_EastBay" = "Red Drum (EB)",
                  "BlueCrabSmall_EastBay" = "Blue Crab  (EB)",
                  "SpottedSeatrout_EastBay" = "Spotted Seatrout (EB)",
                  "Atlanticcroaker_EastBay" = "Atlantic Croaker (EB)",
                  "RedDrum_TrinityBay" = "Red Drum (TB)",
                  "BlueCrabSmall_TrinityBay" = "Blue Crab  (TB)",
                  "SpottedSeatrout_TrinityBay" = "Spotted Seatrout (TB)",
                  "Atlanticcroaker_TrinityBay" = "Atlantic Croaker (TB)",
                  "RedDrum_GalvestonBay" = "Red Drum (GB)",
                  "BlueCrabSmall_GalvestonBay" = "Blue Crab (GB)",
                  "SpottedSeatrout_GalvestonBay" = "Spotted Seatrout (GB)",
                  "Atlanticcroaker_GalvestonBay" = "Atlantic Croaker (GB)"
    )
  )

loodf_AB_Sciaenid_PDSI_forplotting_filtered  <- loodf_AB_Sciaenid_PDSI_forplotting%>%
  filter(!Var2 %in% c("PDSI", "Salinity_CopanoBay", "Salinity_MesquiteBay", "Salinity_AransasBay")) %>%
  mutate(
    Var2 = recode(Var2,
                  "RedDrum_CopanoBay" = "Red Drum (CB)",
                  "BlueCrabSmall_CopanoBay" = "Blue Crab  (CB)",
                  "SpottedSeatrout_CopanoBay" = "Spotted Seatrout (CB)",
                  "Atlanticcroaker_CopanoBay" = "Atlantic Croaker (CB)",
                  "RedDrum_MesquiteBay" = "Red Drum (MB)",
                  "BlueCrabSmall_MesquiteBay" = "Blue Crab  (MB)",
                  "SpottedSeatrout_MesquiteBay" = "Spotted Seatrout (MB)",
                  "Atlanticcroaker_MesquiteBay" = "Atlantic Croaker (MB)",
                  "RedDrum_AransasBay" = "Red Drum (AB)",
                  "BlueCrabSmall_AransasBay" = "Blue Crab  (AB)",
                  "SpottedSeatrout_AransasBay" = "Spotted Seatrout (AB)",
                  "Atlanticcroaker_AransasBay" = "Atlantic Croaker (AB)",
    )
  )

# Create the sciaenid bar plots
order_AB_sciaenid <- c(
  "Atlantic Croaker (AB)",
  "Atlantic Croaker (MB)",
  "Atlantic Croaker (CB)",
  "Blue Crab  (AB)",
  "Blue Crab  (MB)",
  "Blue Crab  (CB)",
  "Red Drum (AB)",
  "Red Drum (MB)",
  "Red Drum (CB)",
  "Spotted Seatrout (AB)",
  "Spotted Seatrout (MB)",
  "Spotted Seatrout (CB)")

PDSI_plot_AB_sciaenid<- ggplot(loodf_AB_Sciaenid_PDSI_forplotting_filtered, aes(x = Var2, y = mean_obs, fill = weather)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(
    ymin = mean_obs - se_obs,
    ymax = mean_obs + se_obs 
  ), width = 0.25, position = position_dodge(0.7)) +  
    scale_fill_manual(values = c("Severely Dry" = "#e9c46a", "Severely Wet" = "#264653")) +
  theme_bw() +  
  labs(x = NULL,y = "Mean Standardized CPUE", fill = NULL, title = "Mission-Aransas Estuary - Sciaenid System")+
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 6),  
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8.5),
    legend.position="none") +
  scale_x_discrete(limits = order_AB_sciaenid) 
print(PDSI_plot_AB_sciaenid)

order_GB_sciaenid <- c(
  "Atlantic Croaker (WB)",
  "Atlantic Croaker (GB)",
  "Atlantic Croaker (EB)",
  "Atlantic Croaker (TB)",
  "Blue Crab  (WB)",
  "Blue Crab (GB)",
  "Blue Crab  (EB)",
  "Blue Crab  (TB)",
  "Red Drum (WB)",
  "Red Drum (GB)",
  "Red Drum (EB)",
  "Red Drum (TB)",
  "Spotted Seatrout (WB)",
  "Spotted Seatrout (GB)",
  "Spotted Seatrout (EB)",
  "Spotted Seatrout (TB)")

PDSI_plot_GB_sciaenid <- ggplot(loodf_GB_Sciaenid_PDSI_forplotting_filtered, aes(x = Var2, y = mean_obs, fill = weather)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(
    ymin = mean_obs - se_obs,
    ymax = mean_obs + se_obs 
  ), width = 0.25, position = position_dodge(0.7)) +  
  scale_fill_manual(values = c("Severely Dry" = "#e9c46a", "Severely Wet" = "#264653")) +
  theme_bw() +  
  labs(x = NULL,y = "Mean Standardized CPUE", fill = NULL, title = "Trinity-San Jacinto Estuary - Sciaenid System")+
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 6),  
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8.5),
    legend.position="none") +
  scale_x_discrete(limits = order_GB_sciaenid) 
print(PDSI_plot_GB_sciaenid)


# Create the keystone bar plots

order_AB_keystone <- c(
  "Mullet (AB)",
  "Mullet (MB)",
  "Mullet (CB)",
  "Menhaden (AB)",
  "Menhaden (MB)",
  "Menhaden (CB)",
  "Alligator Gar (AB)",
  "Alligator Gar (MB)",
  "Alligator Gar (CB)",
  "Bull Shark (AB)",
  "Bull Shark (MB)",
  "Bull Shark (CB)")

PDSI_plot_AB_keystone<- ggplot(loodf_AB_Pred_PDSI_forplotting_filtered, aes(x = Var2, y = mean_obs, fill = weather)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(
    ymin = mean_obs - se_obs,
    ymax = mean_obs + se_obs 
  ), width = 0.25, position = position_dodge(0.7)) +  
  scale_fill_manual(values = c("Severely Dry" = "#e9c46a", "Severely Wet" = "#264653")) +
  theme_bw() +  
  labs(x = NULL,y = "Mean Standardized CPUE", fill = NULL, title = "Mission-Aransas Estuary - Keystone Predator System")+
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 6),  
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8.5),
    legend.position="none")  +
  annotate("text", 
           x = 4.8,  
           y = 1.75, 
           label = paste("MAE Severely Dry Mean PDSI = ", round(mean_PDSI_dry_AB, 1)),
           size = 1.8, hjust = 0, vjust = 0.5, color = "black") +
  annotate("text", 
           x = 4.8,  
           y = 1.55 ,
           label = paste("MAE Severely Wet Mean PDSI = ", round(mean_PDSI_wet_AB, 1)),
           size = 1.8, hjust = 0, vjust = 0.5, color = "black") +
  scale_x_discrete(limits = order_AB_keystone) 
print(PDSI_plot_AB_keystone)

order_GB_keystone <- c(
  "Mullet (WB)",
  "Mullet (GB)",
  "Mullet (EB)",
  "Mullet (TB)",
  "Menhaden (WB)",
  "Menhaden (GB)",
  "Menhaden (EB)",
  "Menhaden (TB)",
  "Alligator Gar (WB)",
  "Alligator Gar (GB)",
  "Alligator Gar (EB)",
  "Alligator Gar (TB)",
  "Bull Shark (WB)",
  "Bull Shark (GB)",
  "Bull Shark (EB)",
  "Bull Shark (TB)")

PDSI_plot_GB_keystone <- ggplot(loodf_GB_Pred_PDSI_forplotting_filtered, aes(x = Var2, y = mean_obs, fill = weather)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(
    ymin = mean_obs - se_obs,
    ymax = mean_obs + se_obs 
  ), width = 0.25, position = position_dodge(0.7)) +  
  scale_fill_manual(values = c("Severely Dry" = "#e9c46a", "Severely Wet" = "#264653")) +
  theme_bw() +  
  labs(x = NULL,y = "Mean Standardized CPUE", fill = NULL, title = "Trinity-San Jacinto Estuary - Keystone Predator System")+
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 6),  
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8.5),
    legend.position = c(0.21, 0.85),  
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key.size = unit(0.3, "cm") 
    ) +
  annotate("text", 
           x = 7.0,  
           y = 2.75,  
           label = paste("TSJE Severely Dry Mean PDSI = ", round(mean_PDSI_dry_GB, 1)),
           size = 1.8, hjust = 0, vjust = 0.5, color = "black") +
  annotate("text", 
           x = 7.0,  
           y = 2.45,  
           label = paste("TSJE Severely Wet Mean PDSI = ", round(mean_PDSI_wet_GB, 1)),
           size = 1.8, hjust = 0, vjust = 0.5, color = "black") +
  scale_x_discrete(limits = order_GB_keystone) 
print(PDSI_plot_GB_keystone)

# export plots
PDSI_Plots_All <- grid.arrange(PDSI_plot_AB_sciaenid, PDSI_plot_AB_keystone, PDSI_plot_GB_sciaenid, PDSI_plot_GB_keystone,
                               ncol = 2, nrow = 2)

ggsave("PDSI_Plots_All.tiff", PDSI_Plots_All, dpi = 300, bg = "white",
       width = 170,
       height = 150,
       units = "mm") 

PDSI_Plots_AB <- grid.arrange(PDSI_plot_AB_sciaenid, PDSI_plot_AB_keystone,
                                    ncol = 1, nrow = 2)

ggsave("PDSI_Plots_AB.png", PDSI_Plots_AB, dpi = 150, bg = "white",
       width = 1600,
       height = 2000,
       units = "px") 

PDSI_Plots_GB <- grid.arrange(PDSI_plot_GB_sciaenid, PDSI_plot_GB_keystone,
                                         ncol = 1, nrow = 2)

ggsave("PDSI_Plots_GB.png", PDSI_Plots_GB, dpi = 300, bg = "white",
       width = 170,
       height = 170,
       units = "mm") 



