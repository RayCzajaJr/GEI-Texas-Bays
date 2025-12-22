

library(dplyr)
library(tidyr)
library(readxl)  
library(gridExtra)
library(ggplot2)
library(readxl)
library(httr)
library(jsonlite)
library(readxl)
library(httr)
library(jsonlite)
library(git2r)
library(tidyverse)
library(patchwork)
library(cowplot)
library(ggspatial)


########## Process (get means per space and time, log transform, z score transform, etc) all raw data 
########## Make basic maps showing spatial dilineations
########## Make time series plots of raw annual values for each species, trophic system and bay



# Get all files from github
repo_api_url <- "https://api.github.com/repos/RayCzajaJr/GEI-Texas-Bays/contents/Input"

response <- GET(repo_api_url)
file_data <- content(response, "text", encoding = "UTF-8")
file_list <- fromJSON(file_data)

excel_files <- file_list$name[grep("\\.xlsx$", file_list$name)]

github_base <- "https://raw.githubusercontent.com/RayCzajaJr/GEI-Texas-Bays/main/Input/"

for (file in excel_files) {
  file_url <- paste0(github_base, file)  
  
  temp_file <- tempfile(fileext = ".xlsx")  
  GET(file_url, write_disk(temp_file, overwrite = TRUE))  
  
  df_name <- gsub("\\.xlsx$", "", file)  
  assign(df_name, read_excel(temp_file), envir = .GlobalEnv)  
}

print("All Excel files have been downloaded and loaded into R.")



# Or get all files locally
setwd("~/Desktop/Winter 2025 DSEM")  
file_list <- list.files("Data", pattern = "\\.xlsx$", full.names = TRUE)
for (file in file_list) {
  df_name <- gsub("\\.xlsx$", "", basename(file))  
  assign(df_name, read_excel(file), envir = .GlobalEnv)
}



# process data for aransas bay sciaenid trophic system first

# merging temp and salinity data from three surveys (Probaly dont need temp data but lets keep it here just in case)
AransasBay_Abiotics <- bind_rows(
  AransasBay_BS_shrimp %>%
    select(MINOR_AREA_CODE, STATION_CODE, START_TEMPERATURE_NUM, START_SALINITY_NUM, DATE),
  AransasBay_GN_cpue %>%
    select(MINOR_AREA_CODE, STATION_CODE, START_TEMPERATURE_NUM, START_SALINITY_NUM, DATE),
  AransasBay_BT %>%
    select(MINOR_AREA_CODE, STATION_CODE, START_TEMPERATURE_NUM, START_SALINITY_NUM, DATE)
)

# function to get data in mean annual values per minor bay in wide format for dsem
processforDSEM <- function(data, species, minor_area_mapping, date_column, result_df) {
  # Step 1: Assign values to the new column MINOR_BAY based on MINOR_AREA_CODE
  data <- data %>%
    mutate(MINOR_BAY = case_when(
      MINOR_AREA_CODE %in% minor_area_mapping$Copano_Bay ~ "Copano Bay",
      MINOR_AREA_CODE %in% minor_area_mapping$Aransas_Bay ~ "Aransas Bay",
      MINOR_AREA_CODE %in% minor_area_mapping$Mesquite_Bay ~ "Mesquite Bay",
      TRUE ~ NA_character_
    ))
  
  # Step 2: Extract YEAR from the specified date column
  data$YEAR <- substr(data[[date_column]], 1, 4)
  
  # Step 3: Calculate the mean value of each species per YEAR and MINOR_BAY
  wide_data <- data %>%
    group_by(YEAR, MINOR_BAY) %>%
    summarize(across(all_of(species), ~ mean(.x, na.rm = TRUE), .names = "{.col}"), .groups = "drop") %>%
    pivot_wider(
      names_from = MINOR_BAY,
      values_from = all_of(species),
      names_glue = "{.value}_{MINOR_BAY}"  # Ensure consistent naming: species_minor_bay
    )
  
  # Step 4: Merge the output with the existing result dataframe
  result_df <- result_df %>%
    full_join(wide_data, by = "YEAR")
  
  return(result_df)
}

# define mapping of MINOR_AREA_CODEs to MINOR_BAY
minor_area_mapping <- list(
  Copano_Bay = c(270, 317, 293, 120, 240, 310),
  Aransas_Bay = c(227, 20, 13, 280, 285, 95, 94),
  Mesquite_Bay = c(143, 70, 315, 250, 90)
)

# make initial df
AransasBay_Sciaenid_Wide <- data.frame(YEAR = character())

# use function to get values for fish species from gill net surveys
AransasBay_Sciaenid_Wide <- processforDSEM(
  data = AransasBay_GN_cpue,
  species = c("RedDrum", "SpottedSeatrout"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Sciaenid_Wide
)


# use function to get values for juvenile blue crabs from bag seine surveys
AransasBay_Sciaenid_Wide <- processforDSEM(
  data = AransasBay_BS_bluecrabs,
  species = c("BlueCrabSmall"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Sciaenid_Wide
)

# use function to get values for croaker from bag seine surveys
AransasBay_Sciaenid_Wide <- processforDSEM(
  data = AransasBay_BS_Croaker,
  species = c("Atlantic croaker"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Sciaenid_Wide
)

# use function to get values for salinity from all gill net, bag seine and bottom trawl surveys
AransasBay_Sciaenid_Wide <- processforDSEM(
  data = AransasBay_Abiotics,
  species = c("START_SALINITY_NUM"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Sciaenid_Wide
)


# the df now has all the needed data (I think?). we still need to transform/standardize the values, but first lets clean up the df (remove and rename some columns)
AransasBay_Sciaenid_Wide <- AransasBay_Sciaenid_Wide %>%
  select(-START_SALINITY_NUM_NA, -BlueCrabSmall_NA, -"Atlantic croaker_NA")

AransasBay_Sciaenid_Wide <- AransasBay_Sciaenid_Wide %>%
  rename_with(
    .fn = ~ gsub("^START_SALINITY_NUM_", "Salinity_", .),
    .cols = starts_with("START_SALINITY_NUM_")
  )

# remove spaces from column names
colnames(AransasBay_Sciaenid_Wide) <- gsub(" ", "", colnames(AransasBay_Sciaenid_Wide))


# okay now lets transform the values and create a new df
AransasBay_Sciaenid_Wide_Trans <- AransasBay_Sciaenid_Wide %>%
  # Step 1: Apply log transformation
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ log(. + 1)  # Log transformation with a small shift to avoid log(0)
  )) %>%
  # Step 2: Apply standardization (mean = 0, SD = 1)
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ (. - mean(.)) / sd(.)
  ))

# get PDSI values in there
PDSI_annual_ara$YEAR <- as.character(PDSI_annual_ara$YEAR)
AransasBay_Sciaenid_Wide_Trans <- AransasBay_Sciaenid_Wide_Trans %>%
  left_join(PDSI_annual_ara %>% select(YEAR, PDSI), by = "YEAR")


# this data (sciaenid system - aransas bay) is now ready for dsem

# now lets process more data (keystone predator system - aransas bay)

# make initial df
AransasBay_Pred_Wide <- data.frame(YEAR = character())

# use function to get values for bullshark from gill net surveys
AransasBay_Pred_Wide <- processforDSEM(
  data = AransasBay_GN_BullShark,
  species = c("BullShark"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Pred_Wide
)

# use function to get values for gar from gill net surveys
AransasBay_Pred_Wide <- processforDSEM(
  data = AransasBay_GN_Gar,
  species = c("AlligatorGar"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Pred_Wide
)

# renaming menhaden and mullet columns 
AransasBay_GN_cpue <- AransasBay_GN_cpue %>%
  rename(AllMenhaden = GulfMenhaden)

AransasBay_GN_Mullet <- AransasBay_GN_Mullet %>%
  rename(AllMullet = StripedMullet)


# use function to get values for mullet from gill net surveys
AransasBay_Pred_Wide <- processforDSEM(
  data = AransasBay_GN_Mullet,
  species = c("AllMullet"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Pred_Wide
)

# use function to get values for menhaden from gill net surveys
AransasBay_Pred_Wide <- processforDSEM(
  data = AransasBay_GN_cpue,
  species = c("AllMenhaden"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Pred_Wide
)

# use function to get values for salinity from all gill net, bag seine and bottom trawl surveys
AransasBay_Pred_Wide <- processforDSEM(
  data = AransasBay_Abiotics,
  species = c("START_SALINITY_NUM"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Pred_Wide
)

# the df now has all the needed data (I think?). we still need to transform/standardize the values, but first lets clean up the df (remove and rename some columns)
AransasBay_Pred_Wide <- AransasBay_Pred_Wide %>%
  select(-START_SALINITY_NUM_NA)

AransasBay_Pred_Wide <- AransasBay_Pred_Wide %>%
  rename_with(
    .fn = ~ gsub("^START_SALINITY_NUM_", "Salinity_", .),
    .cols = starts_with("START_SALINITY_NUM_")
  )

# quick inspection of some temporal trends
time_series_data <- AransasBay_Pred_Wide %>%
  select(YEAR, starts_with("BullShark")) %>%
  pivot_longer(
    cols = -YEAR,
    names_to = c("Species", "Bay"),
    names_sep = "_",  # Assumes the column names are in "Species_Bay" format
    values_to = "Mean_CPUE"
  )

ggplot(time_series_data, aes(x = YEAR, y = Mean_CPUE, color = Bay, group = interaction(Species, Bay))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Species, scales = "free_y") +
  labs(
    title = "Time Series of BullShark",
    x = "Year",
    y = "Mean CPUE",
    color = "Bay"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# remove spaces from column names
colnames(AransasBay_Pred_Wide) <- gsub(" ", "", colnames(AransasBay_Pred_Wide))


# okay now lets transform the values and create a new df
AransasBay_Pred_Wide_Trans <- AransasBay_Pred_Wide %>%
  # Step 1: Apply log transformation
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ log(. + 1)  # Log transformation with a small shift to avoid log(0)
  )) %>%
  # Step 2: Apply standardization (mean = 0, SD = 1)
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ (. - mean(.)) / sd(.)
  ))

# get PDSI values in there
PDSI_annual_ara$YEAR <- as.character(PDSI_annual_ara$YEAR)
AransasBay_Pred_Wide_Trans <- AransasBay_Pred_Wide_Trans %>%
  left_join(PDSI_annual_ara %>% select(YEAR, PDSI), by = "YEAR")




#################################################################################

# now galveston bay and the scieanid trophic system

# merging temp and salinity data from three surveys (Probaly dont need temp data but lets keep it here just in case)
GalvestonBay_Abiotics <- bind_rows(
  GalvestonBay_BS_shrimp %>%
    select(MINOR_AREA_CODE, STATION_CODE, START_TEMPERATURE_NUM, START_SALINITY_NUM, DATE),
  GalvestonBay_GN_cpue %>%
    select(MINOR_AREA_CODE, STATION_CODE, START_TEMPERATURE_NUM, START_SALINITY_NUM, DATE),
  GalvestonBay_BT_cpue %>%
    select(MINOR_AREA_CODE, STATION_CODE, START_TEMPERATURE_NUM, START_SALINITY_NUM, DATE)
)

# function to get data in mean annual values per minor bay in wide format for dsem
processforDSEM_GB <- function(data, species, minor_area_mapping, date_column, result_df) {
  # Step 1: Assign values to the new column MINOR_BAY based on MINOR_AREA_CODE
  data <- data %>%
    mutate(MINOR_BAY = case_when(
      MINOR_AREA_CODE %in% minor_area_mapping$Galveston_Bay ~ "Galveston Bay",
      MINOR_AREA_CODE %in% minor_area_mapping$Trinity_Bay ~ "Trinity Bay",
      MINOR_AREA_CODE %in% minor_area_mapping$West_Bay ~ "West Bay",
      MINOR_AREA_CODE %in% minor_area_mapping$East_Bay ~ "East Bay",
      TRUE ~ NA_character_
    ))
  
  # Step 2: Extract YEAR from the specified date column
  data$YEAR <- substr(data[[date_column]], 1, 4)
  
  # Step 3: Calculate the mean value of each species per YEAR and MINOR_BAY
  wide_data <- data %>%
    group_by(YEAR, MINOR_BAY) %>%
    summarize(across(all_of(species), ~ mean(.x, na.rm = TRUE), .names = "{.col}"), .groups = "drop") %>%
    pivot_wider(
      names_from = MINOR_BAY,
      values_from = all_of(species),
      names_glue = "{.value}_{MINOR_BAY}"  # Ensure consistent naming: species_minor_bay
    )
  
  # Step 4: Merge the output with the existing result dataframe
  result_df <- result_df %>%
    full_join(wide_data, by = "YEAR")
  
  return(result_df)
}

# define mapping of MINOR_AREA_CODEs to MINOR_BAY
minor_area_mapping_GB <- list(
  Galveston_Bay = c(180, 141, 241, 142),
  Trinity_Bay = c(53, 64, 319, 318, 54, 312, 12, 330),
  East_Bay = c(150),
  West_Bay = c(311, 201, 181, 61, 350, 191, 100, 225, 261, 50, 222, 131, 110, 144)
)

# make initial df
GalvestonBay_Sciaenid_Wide <- data.frame(YEAR = character())

# use function to get values for fish species from gill net surveys
GalvestonBay_Sciaenid_Wide <- processforDSEM_GB(
  data = GalvestonBay_GN_cpue,
  species = c("RedDrum", "SpottedSeatrout"),
  minor_area_mapping = minor_area_mapping_GB,
  date_column = "DATE",
  result_df = GalvestonBay_Sciaenid_Wide
)

# use function to get values for croaker from bag seine surveys
GalvestonBay_Sciaenid_Wide <- processforDSEM_GB(
  data = GalvestonBay_BS_Croaker,
  species = c("Atlantic croaker"),
  minor_area_mapping = minor_area_mapping_GB,
  date_column = "DATE",
  result_df = GalvestonBay_Sciaenid_Wide
)


# use function to get values for juvenile blue crabs from bag seine surveys
GalvestonBay_Sciaenid_Wide <- processforDSEM_GB(
  data = GalvestonBay_BS_bluecrabs,
  species = c("BlueCrabSmall"),
  minor_area_mapping = minor_area_mapping_GB,
  date_column = "DATE",
  result_df = GalvestonBay_Sciaenid_Wide
)

# use function to get values forsalinity from all gill net, bag seine and bottom trawl surveys
GalvestonBay_Sciaenid_Wide <- processforDSEM_GB(
  data = GalvestonBay_Abiotics,
  species = c("START_SALINITY_NUM"),
  minor_area_mapping = minor_area_mapping_GB,
  date_column = "DATE",
  result_df = GalvestonBay_Sciaenid_Wide
)


# the df now has all the needed data (I think?). we still need to transform/standardize the values, but first lets clean up the df (remove and rename some columns)
GalvestonBay_Sciaenid_Wide <- GalvestonBay_Sciaenid_Wide %>%
  select(-START_SALINITY_NUM_NA, -BlueCrabSmall_NA, -"Atlantic croaker_NA")

GalvestonBay_Sciaenid_Wide <- GalvestonBay_Sciaenid_Wide %>%
  rename_with(
    .fn = ~ gsub("^START_SALINITY_NUM_", "Salinity_", .),
    .cols = starts_with("START_SALINITY_NUM_")
  )


# remove spaces from column names
colnames(GalvestonBay_Sciaenid_Wide) <- gsub(" ", "", colnames(GalvestonBay_Sciaenid_Wide))


# okay now lets transform the values and create a new df
GalvestonBay_Sciaenid_Wide_Trans <- GalvestonBay_Sciaenid_Wide %>%
  # Step 1: Apply log transformation
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ log(. + 1)  # Log transformation with a small shift to avoid log(0)
  )) %>%
  # Step 2: Apply standardization (mean = 0, SD = 1)
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)  # Standardization
  ))

# get PDSI values in there
PDSI_annual_gal$YEAR <- as.character(PDSI_annual_gal$YEAR)
GalvestonBay_Sciaenid_Wide_Trans <- GalvestonBay_Sciaenid_Wide_Trans %>%
  left_join(PDSI_annual_gal %>% select(YEAR, PDSI), by = "YEAR")



# this data (sciaenid system - Galveston bay) is now ready for dsem

# now lets process more data (keystone predator system - Galveston bay)

# make initial df
GalvestonBay_Pred_Wide <- data.frame(YEAR = character())

# use function to get values for bullshark from gill net surveys
GalvestonBay_Pred_Wide <- processforDSEM_GB(
  data = GalvestonBay_GN_BullShark,
  species = c("BullShark"),
  minor_area_mapping = minor_area_mapping_GB,
  date_column = "DATE",
  result_df = GalvestonBay_Pred_Wide
)

# use function to get values for gar from gill net surveys
GalvestonBay_Pred_Wide <- processforDSEM_GB(
  data = GalvestonBay_GN_Gar,
  species = c("AlligatorGar"),
  minor_area_mapping = minor_area_mapping_GB,
  date_column = "DATE",
  result_df = GalvestonBay_Pred_Wide
)

# renaming menhaden and mullet column names
GalvestonBay_GN_cpue <- GalvestonBay_GN_cpue %>%
  rename(Menhaden = GulfMenhaden)

GalvestonBay_GN_Mullet <- GalvestonBay_GN_Mullet %>%
  rename(Mullet = StripedMullet)


# use function to get values for mullet from gill net surveys
GalvestonBay_Pred_Wide <- processforDSEM_GB(
  data = GalvestonBay_GN_Mullet,
  species = c("Mullet"),
  minor_area_mapping = minor_area_mapping_GB,
  date_column = "DATE",
  result_df = GalvestonBay_Pred_Wide
)

# use function to get values for mullet from gill net surveys
GalvestonBay_Pred_Wide <- processforDSEM_GB(
  data = GalvestonBay_GN_cpue,
  species = c("Menhaden"),
  minor_area_mapping = minor_area_mapping_GB,
  date_column = "DATE",
  result_df = GalvestonBay_Pred_Wide
)


# use function to get values for salinity from all gill net, bag seine and bottom trawl surveys
GalvestonBay_Pred_Wide <- processforDSEM_GB(
  data = GalvestonBay_Abiotics,
  species = c("START_SALINITY_NUM"),
  minor_area_mapping = minor_area_mapping_GB,
  date_column = "DATE",
  result_df = GalvestonBay_Pred_Wide
)

# the df now has all the needed data (I think?). we still need to transform/standardize the values, but first lets clean up the df (remove and rename some columns)
GalvestonBay_Pred_Wide <- GalvestonBay_Pred_Wide %>%
  select(-START_SALINITY_NUM_NA)

GalvestonBay_Pred_Wide <- GalvestonBay_Pred_Wide %>%
  rename_with(
    .fn = ~ gsub("^START_SALINITY_NUM_", "Salinity_", .),
    .cols = starts_with("START_SALINITY_NUM_")
  )


# remove spaces from column names
colnames(GalvestonBay_Pred_Wide) <- gsub(" ", "", colnames(GalvestonBay_Pred_Wide))

# okay now lets transform the values and create a new df
GalvestonBay_Pred_Wide_Trans <- GalvestonBay_Pred_Wide %>%
  # okay now lets transform the values and create a new df
  # Step 1: Apply log transformation
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ log(. + 1)  # Log transformation with a small shift to avoid log(0)
  )) %>%
  # Step 2: Apply standardization (mean = 0, SD = 1)
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)  # Standardization
  ))

# get PDSI values in there
PDSI_annual_gal$YEAR <- as.character(PDSI_annual_gal$YEAR)
GalvestonBay_Pred_Wide_Trans <- GalvestonBay_Pred_Wide_Trans %>%
  left_join(PDSI_annual_gal %>% select(YEAR, PDSI), by = "YEAR")




###### MAP TIME 



# make maps out of minor bays and clusters
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# extract base data for map
coastline_10 <- ne_coastline(scale = 10, returnclass = "sf")
world <- ne_countries(scale = 10, returnclass = "sf")
us_states <- ne_states(country = "United States of America", returnclass = "sf")  

### GB

# define mapping of MINOR_AREA_CODEs to MINOR_BAY
minor_area_mapping_GB <- list(
  Galveston_Bay = c(180, 141, 241, 142),
  Trinity_Bay = c(53, 64, 319, 318, 54, 312, 12, 330),
  East_Bay = c(150),
  West_Bay = c(311, 201, 181, 61, 350, 191, 100, 225, 261, 50, 222, 131, 110, 144, 150)
)

# add new column "Bay_Area" to classify each row based on MINOR_AREA_CODE
GalvestonBay_GN_cpue_map <- GalvestonBay_GN_cpue %>%
  mutate(Bay_Area = case_when(
    MINOR_AREA_CODE %in% minor_area_mapping_GB$Galveston_Bay ~ "Galveston Bay",
    MINOR_AREA_CODE %in% minor_area_mapping_GB$Trinity_Bay ~ "Trinity Bay",
    MINOR_AREA_CODE %in% minor_area_mapping_GB$East_Bay ~ "East Bay",
    MINOR_AREA_CODE %in% minor_area_mapping_GB$West_Bay ~ "West Bay",
    TRUE ~ "Other"  # Optional: Handles any unexpected codes
  ))

# define color palette
bay_colors_GB <- c("Galveston Bay" = "#e9c46a", 
                "Trinity Bay" = "#f4a261", 
                "East Bay" = "peachpuff3",
                "West Bay" = "#f28482")

GBMAP<-ggplot() +
  geom_sf(data = world, color = "black") +
  coord_sf(xlim = c(-95.4, -94.5),
           ylim = c(29.0, 29.8)) +
  geom_point(data = GalvestonBay_GN_cpue_map, aes(x = X, y = Y, color = Bay_Area), size = 1.5) +
  labs(x = "Longitude", y = "Latitude", title = "Trinity-San Jacinto Estuary", color = "Bay") +
  theme_bw()+
  scale_color_manual(values = bay_colors_GB) +
  theme(
    legend.position = "bottom",  
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 8), 
    axis.text = element_text(size = 6), 
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),    # reduce margin around legend
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0), # reduce margin for the box
    legend.key.height = unit(0.3, "cm"),                   # shrink keys
    legend.key.width = unit(0.6, "cm"),
    axis.title = element_text(size = 6),  
    plot.title = element_text(size = 10)  
  )+
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) + 
  annotation_scale(location = "tl", width_hint = 0.3, text_cex = 0.7, unit_category = "metric") + 
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.9, "in"),
                         style = north_arrow_fancy_orienteering)

# AB 

# define mapping of MINOR_AREA_CODEs to MINOR_BAY
minor_area_mapping_AB <- list(
  Copano_Bay = c(270, 317, 293, 120, 240, 310),
  Aransas_Bay = c(227, 20, 13, 280, 285, 95, 94),
  Mesquite_Bay = c(143, 70, 315, 250, 90)
)

# add new column "Bay_Area" to classify each row based on MINOR_AREA_CODE
AransasBay_GN_cpue_map <- AransasBay_GN_cpue %>%
  mutate(Bay_Area = case_when(
    MINOR_AREA_CODE %in% minor_area_mapping_AB$Copano_Bay ~ "Copano Bay",
    MINOR_AREA_CODE %in% minor_area_mapping_AB$Aransas_Bay ~ "Aransas Bay",
    MINOR_AREA_CODE %in% minor_area_mapping_AB$Mesquite_Bay ~ "Mesquite Bay",
    TRUE ~ "Other"  # Optional: Handles any unexpected codes
  ))

# define color palette
bay_colors_AB <- c("Copano Bay" = "#264653", 
                   "Aransas Bay" =  "#2a9d8f", 
                   "Mesquite Bay" = "#90a955")

ABMAP<-ggplot() +
  geom_sf(data = world, color = "black") +
  coord_sf(xlim = c(-97.3, -96.7),
           ylim = c(27.8, 28.3)) +
  geom_point(data = AransasBay_GN_cpue_map, aes(x = X, y = Y, color = Bay_Area), size = 1.5) +
  labs(x = "Longitude", y = "Latitude", title = "Mission-Aransas Estuary", color = "Bay") +
  theme_bw()+
  scale_color_manual(values = bay_colors_AB) +
  theme(
    legend.position = "bottom",  
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 8), 
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),    # reduce margin around legend
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0), # reduce margin for the box
    legend.key.height = unit(0.3, "cm"),                   # shrink keys
    legend.key.width = unit(0.6, "cm"),
    axis.text = element_text(size = 6), 
    axis.title = element_text(size = 6),  
    plot.title = element_text(size = 10)  
  )+
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) + 
  annotation_scale(location = "tl", width_hint = 0.3, text_cex = 0.7, unit_category = "metric") 
  

# make inset for whole state of texas
TEXMAP<-ggplot() +
  geom_sf(data = world, color = "black") +
  geom_sf(data = us_states, color = "black", fill = NA, linetype = "dashed") +  
  coord_sf(xlim = c(-106.3, -93.7),
           ylim = c(25.8, 36.3)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()+
  theme(
    axis.text = element_blank(),  
    axis.title = element_blank(),  
    axis.ticks = element_blank()   
  )

TEXMAP_grob <- ggplotGrob(TEXMAP)

GBMAP <- GBMAP +
  annotation_custom(grob = TEXMAP_grob, xmin = -94.79, xmax = -94.5, ymin = 28.96, ymax = 29.3)
ABMAP <- ABMAP +
  annotation_custom(grob = TEXMAP_grob, xmin = -96.89, xmax = -96.69, ymin = 27.79, ymax = 27.99)

# combine the three maps into one layout
combined_map <- grid.arrange(ABMAP, GBMAP, ncol = 1, nrow = 2)
ggsave("combined_map.tiff", combined_map, width = 85, height = 190, dpi = 300, units = "mm")




####### Make time series plots of raw annual values for each species, trophic system and bay




# Reusable function for plotting
plot_species_timeseries <- function(df_wide, bay_colors, y_labels) {
  
  # Make sure YEAR is numeric
  df_wide <- df_wide %>%
    mutate(YEAR = as.numeric(YEAR))
  
  # Long format
  long_df <- df_wide %>%
    pivot_longer(cols = -YEAR, names_to = "variable", values_to = "value") %>%
    separate(variable, into = c("species", "bay"), sep = "_")
  
  # Drop rows where species isn't in y_labels
  long_df <- long_df %>%
    filter(species %in% names(y_labels))
  
  # Order species based on y_labels order
  long_df <- long_df %>%
    mutate(species = factor(species, levels = names(y_labels)))
  
  # Custom theme
  custom_theme <- theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 14),
      plot.title = element_blank()
    )
  
  # Create plots
  species_plots <- map(levels(long_df$species), function(sp) {
    df_sp <- filter(long_df, species == sp)
    
    ggplot(df_sp, aes(x = YEAR, y = value, color = bay)) +
      geom_line(size = 1.5) +
      scale_color_manual(values = bay_colors) +
      labs(
        y = y_labels[[sp]],
        x = "Year"
      ) +
      custom_theme
  })
  
  # Create large legend
  legend_plot <- ggplot(long_df, aes(x = YEAR, y = value, color = bay)) +
    geom_line(size = 2) +
    scale_color_manual(values = bay_colors) +
    labs(color = "Bay") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.key.size = unit(1.5, "lines")
    )
  
  legend <- get_legend(legend_plot)
  
  # Combine plots + legend in 3x2 grid
  final_plot <- (
    species_plots[[1]] + species_plots[[2]] + species_plots[[3]] +
      species_plots[[4]] + species_plots[[5]] + wrap_elements(legend)
  ) +
    plot_layout(ncol = 3)
  
  return(final_plot)
}

# ---- Aransas Sciaenid ----
aransas_colors <- c(
  "CopanoBay" = "#264653",
  "AransasBay" = "#2a9d8f",
  "MesquiteBay" = "#90a955"
)

galveston_colors <- c(
  "GalvestonBay" = "#e9c46a",
  "TrinityBay" = "#f4a261",
  "EastBay" = "peachpuff3",
  "WestBay" = "#f28482"
)

sciaenid_labels <- c(
  BlueCrabSmall = "Blue Crab Small CPUE",
  Atlanticcroaker = "Atlantic Croaker CPUE",
  Salinity = "Salinity (ppt)",
  RedDrum = "Red Drum CPUE",
  SpottedSeatrout = "Spotted Seatrout CPUE"
)

pred_labels <- c(
  Mullet = "Mullet CPUE",
  Menhaden = "Menhaden CPUE",
  Salinity = "Salinity (ppt)",
  BullShark = "Bull Shark CPUE",
  AlligatorGar = "Alligator Gar CPUE"
)

pred_labels_AB <- c(
  AllMullet = "Mullet CPUE",
  AllMenhaden = "Menhaden CPUE",
  Salinity = "Salinity (ppt)",
  BullShark = "Bull Shark CPUE",
  AlligatorGar = "Alligator Gar CPUE"
)

TimeSeriesAB_Sci<-plot_species_timeseries(AransasBay_Sciaenid_Wide, aransas_colors, sciaenid_labels)
TimeSeriesAB_Pred<-plot_species_timeseries(AransasBay_Pred_Wide, aransas_colors, pred_labels_AB)
TimeSeriesGB_Sci<-plot_species_timeseries(GalvestonBay_Sciaenid_Wide, galveston_colors, sciaenid_labels)
TimeSeriesGB_Pred<-plot_species_timeseries(GalvestonBay_Pred_Wide, galveston_colors, pred_labels)

ggsave("TimeSeriesAB_Sci.png", TimeSeriesAB_Sci, width = 16, height = 9, dpi = 200)
ggsave("TimeSeriesAB_Pred.png", TimeSeriesAB_Pred, width = 16, height = 9, dpi = 200)
ggsave("TimeSeriesGB_Sci.png", TimeSeriesGB_Sci, width = 16, height = 9, dpi = 200)
ggsave("TimeSeriesGB_Pred.png", TimeSeriesGB_Pred, width = 16, height = 9, dpi = 200)

# make time series plots of whole bay means of salinity and PDSI for both major bays 

# Calculate row-wise means for salinity
AransasBaySalinity <- rowMeans(
  AransasBay_Sciaenid_Wide[, c("Salinity_AransasBay", "Salinity_CopanoBay", "Salinity_MesquiteBay")],
  na.rm = TRUE)
GalvestonBaySalinity <- rowMeans(
  GalvestonBay_Sciaenid_Wide[, c("Salinity_WestBay", "Salinity_EastBay", "Salinity_GalvestonBay", "Salinity_TrinityBay")],
  na.rm = TRUE)

YEAR <- AransasBay_Sciaenid_Wide_Trans$YEAR
AransasBayPDSI <- AransasBay_Sciaenid_Wide_Trans$PDSI
GalvestonBayPDSI <- GalvestonBay_Sciaenid_Wide_Trans$PDSI

AbioticMeans <- data.frame(
  YEAR = YEAR,
  AransasBaySalinity = AransasBaySalinity,
  GalvestonBaySalinity = GalvestonBaySalinity,
  AransasBayPDSI = AransasBayPDSI,
  GalvestonBayPDSI = GalvestonBayPDSI)

AbioticMeans_long <- AbioticMeans %>%
  pivot_longer(
    cols = c(AransasBaySalinity, GalvestonBaySalinity, AransasBayPDSI, GalvestonBayPDSI),
    names_to = c("Bay", "Variable"),
    names_pattern = "(AransasBay|GalvestonBay)(Salinity|PDSI)",
    values_to = "Value" ) %>%
  pivot_wider(
    names_from = Variable,
    values_from = Value)

SALPLOT<-ggplot(AbioticMeans_long, aes(x = YEAR, y = Salinity, color = Bay, group = Bay)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c(
      "AransasBay" = "#2a9d8f",
      "GalvestonBay" = "#f28482"))+
  labs(
    x = NULL, 
    color = "Bay")+
  theme_bw()+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_blank(),
          plot.title = element_blank())

PDSIPLOT<-ggplot(AbioticMeans_long, aes(x = YEAR, y = PDSI, color = Bay, group = Bay)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c(
      "AransasBay" = "#2a9d8f",
      "GalvestonBay" = "#f28482"),
    labels = c(
      "AransasBay" = "Mission-Aransas Estuary",
      "GalvestonBay" = "Trinity-San Jacinto Estuary"))+
  labs(
    x = "Year",
    y = "PDSI",
    color = "Bay")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    legend.position = c(0.05, 0.05),      
    legend.justification = c(0, 0),  
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11),
    legend.background = element_rect(fill = "white", color = NA), 
    plot.title = element_blank())


# combine the three maps into one layout
ABIOTICPLOT <- grid.arrange(SALPLOT, PDSIPLOT, ncol = 1, nrow = 2)
ggsave("ABIOTICPLOT.png", ABIOTICPLOT, width = 7, height = 9, dpi = 200)

mean(AbioticMeans_long$Salinity[AbioticMeans_long$Bay == "AransasBay"], na.rm = TRUE)
mean(AbioticMeans_long$Salinity[AbioticMeans_long$Bay == "GalvestonBay"], na.rm = TRUE)
mean(AbioticMeans_long$PDSI[AbioticMeans_long$Bay == "AransasBay"], na.rm = TRUE)
mean(AbioticMeans_long$PDSI[AbioticMeans_long$Bay == "GalvestonBay"], na.rm = TRUE)



