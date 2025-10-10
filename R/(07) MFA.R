library("FactoMineR")
library("factoextra")

# add new columns in dfs
loodf_AB_Pred_PDSI_forplotting$MajorBay <- "AransasBay"
loodf_AB_Pred_PDSI_forplotting$TrophicSystem <- "Keystone Predator"

loodf_GB_Pred_PDSI_forplotting$MajorBay <- "GalvestonBay"
loodf_GB_Pred_PDSI_forplotting$TrophicSystem <- "Keystone Predator"

loodf_AB_Sciaenid_PDSI_forplotting$MajorBay <- "AransasBay"
loodf_AB_Sciaenid_PDSI_forplotting$TrophicSystem <- "Sciaenid"

loodf_GB_Sciaenid_PDSI_forplotting$MajorBay <- "GalvestonBay"
loodf_GB_Sciaenid_PDSI_forplotting$TrophicSystem <- "Sciaenid"

# combine dfs
loodf_PDSI_forPCA <- rbind(
  loodf_AB_Pred_PDSI_forplotting,
  loodf_GB_Pred_PDSI_forplotting,
  loodf_AB_Sciaenid_PDSI_forplotting,
  loodf_GB_Sciaenid_PDSI_forplotting
)

# clean up new df 
loodf_PDSI_forPCA <- loodf_PDSI_forPCA %>%
  # Remove observations where Var2 starts with "Salinity" or "PDSI"
  filter(!str_starts(Var2, "Salinity"), !str_starts(Var2, "PDSI")) %>%  # Separate Var2 into Species and MinorBay using the underscore
  separate(Var2, into = c("Species", "MinorBay"), sep = "_", remove = FALSE) %>%
  # Replace specific terms in the Species column
  mutate(Species = case_when(
    Species == "BlueCrabSmall" ~ "BlueCrab",
    Species == "AllMullet" ~ "Mullet",
    Species == "AllMenhaden" ~ "Menhaden",
    TRUE ~ Species))


# rearrange the df for mfa
loodf_PDSI_forPCA <- loodf_PDSI_forPCA %>%
  select(
    mean_obs,
    MajorBay,
    TrophicSystem,
    weather
  )

loodf_PDSI_forPCA <- loodf_PDSI_forPCA %>%
  mutate(MajorBay = ifelse(MajorBay == "AransasBay", "Aransas Bay", MajorBay)) %>%
  mutate(MajorBay = ifelse(MajorBay == "GalvestonBay", "Galveston Bay", MajorBay))
  
res.mfa <- MFA(loodf_PDSI_forPCA,
               group = c(1, 1, 1, 1),  
               type = c("s", "n", "n", "n"),
               name.group = c("est", "MajorBay",
                              "TrophicSystem", "weather"),
               graph = FALSE)
res.mfa$group
fviz_contrib(res.mfa, "group", axes = 1)
fviz_contrib(res.mfa, "group", axes = 2)

my_colors_weather <- c(
  "Severely Dry" = "#e9c46a",
  "Severely Wet" = "#264653")

weatherbiplot<-fviz_mfa_ind(res.mfa,
             habillage = "weather",  
             palette = my_colors_weather,
             addEllipses = TRUE,
             mean.point = FALSE,
             label = "none",
             repel = TRUE,
             pointsize = 3,
             legend.title = "Weather Condition")+
  theme(
    legend.position = "right",         # Move the legend to the top
    legend.title = element_text(size = 12),  # Customize legend title font size
    plot.title = element_blank()     # Remove the plot title
  )

my_colors_bay <- c(
  "Aransas Bay" =  "#2a9d8f",
  "Galveston Bay" = "#f28482")

majorbaybiplot<-fviz_mfa_ind(res.mfa,
             habillage = "MajorBay",  
             palette = my_colors_bay,
             addEllipses = TRUE,
             mean.point = FALSE,
             label = "none",
             repel = TRUE,
             pointsize = 3,
             legend.title = "Major Bay")+
  theme(
    legend.position = "right",         # Move the legend to the top
    legend.title = element_text(size = 12),  # Customize legend title font size
    plot.title = element_blank()     # Remove the plot title
  )

biplots <- grid.arrange(majorbaybiplot, weatherbiplot, ncol = 1, nrow = 2)

ggsave("biplots.png", biplots, dpi = 250, bg = "white",
       width = 2000,
       height = 2600,
       units = "px")

###############################
