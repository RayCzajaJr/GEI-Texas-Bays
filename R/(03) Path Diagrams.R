library(dsem)
library(ggplot2)
library(ggpubr)
library(ggraph)
library(phylopath)
library(dplyr)
library(ggdag)
library(readxl)
library(qgraph)
library(DHARMa)


####### Make custom path diagrams for each (lowest AIC) trophic model and major bay 
####### Make 'empty' path diagrams as conceptual models 
####### Make dot plot / heat map with significant paths from all lowest AIC models 


# Define custom layout with (x, y) coordinates for path diagrams
custom_layout <- matrix(c(
  -1,  1.9, # salinity
  1,  1,  # spottedseatrout
  -1,  1.1,   # PDSI
  0, 1.6,   # croaker
  1,  2,   # drum 
  0, 1.4    # blue crab
), ncol = 2, byrow = TRUE)

custom_layout_padded <- custom_layout * 1.5  # 20% bigger

# Function to extract coefficients per bay
get_part_Sciaenid_site <- function(x, site_suffix) {
  vars <- c("Salinity","SpottedSeatrout","PDSI",
            "Atlanticcroaker","RedDrum","BlueCrabSmall")
  index <- sapply(vars, function(v) {
    if (v == "PDSI") {
      m <- grep("^PDSI$", rownames(x$coef))
    } else {
      m <- grep(paste0("^", v, "_", site_suffix, "$"), rownames(x$coef))}
    if (length(m) == 0) return(NA_integer_) else return(m)})
  keep <- !is.na(index)
  index <- index[keep]
  vars  <- vars[keep]
  if (length(index) < 2) return(NULL)
  x$coef <- x$coef[index, index, drop = FALSE]
  dimnames(x$coef) <- list(vars, vars)
  return(x) }

get_part_Pred_site <- function(x, site_suffix) {
  vars <- c("Salinity","AlligatorGar","PDSI",
            "AllMullet","BullShark","AllMenhaden")
  index <- sapply(vars, function(v) {
    if (v == "PDSI") {
      m <- grep("^PDSI$", rownames(x$coef))
    } else {
      m <- grep(paste0("^", v, "_", site_suffix, "$"), rownames(x$coef))}
    if (length(m) == 0) return(NA_integer_) else return(m)})
  keep <- !is.na(index)
  index <- index[keep]
  vars  <- vars[keep]
  if (length(index) < 2) return(NULL)
  x$coef <- x$coef[index, index, drop = FALSE]
  dimnames(x$coef) <- list(vars, vars)
  return(x) }

get_part_Pred_site_GB <- function(x, site_suffix) {
  vars <- c("Salinity","AlligatorGar","PDSI",
            "Mullet","BullShark","Menhaden")
  index <- sapply(vars, function(v) {
    if (v == "PDSI") {
      m <- grep("^PDSI$", rownames(x$coef))
    } else {
      m <- grep(paste0("^", v, "_", site_suffix, "$"), rownames(x$coef))}
    if (length(m) == 0) return(NA_integer_) else return(m)})
  keep <- !is.na(index)
  index <- index[keep]
  vars  <- vars[keep]
  if (length(index) < 2) return(NULL)
  x$coef <- x$coef[index, index, drop = FALSE]
  dimnames(x$coef) <- list(vars, vars)
  return(x) }

# Function to create a path diagram with lags
plot_qgraph_combined_with_lag <- function(data_ts, coef_matrix_yes_lag, plot_title = "Path Diagram", bg_color = "transparent") {
  
  # Extract variable names
  abbrev_names <- colnames(data_ts)
  
  # Create an edge weight matrix initialized with zeros for combined plot
  combined_coef_matrix <- coef_matrix_yes_lag
  combined_coef_matrix[,] <- 0
  
  # Store line types: default solid (1), change to dashed (2) for lagged effects
  line_types <- matrix(1, nrow = nrow(combined_coef_matrix), ncol = ncol(combined_coef_matrix))
  
  # Store curvature settings
  edge_curvature <- matrix(0, nrow = nrow(combined_coef_matrix), ncol = ncol(combined_coef_matrix))
  
  # Populate matrix with lagged effects (dashed lines)
  nonzero_yes_lag <- coef_matrix_yes_lag != 0
  combined_coef_matrix[nonzero_yes_lag] <- coef_matrix_yes_lag[nonzero_yes_lag]
  
  # Apply inward curvature for lagged edges (opposite direction)
  edge_curvature[nonzero_yes_lag] <- -1.5  # Inward curve
  
  # Set line types: solid (1) for non-lagged, dashed (2) for lagged
  line_types[nonzero_yes_lag] <- 2  # Dashed line for lagged effects
  
  # Set plot margins
  par(mar = c(0, 13, 0, 0))  
  
  # Plot the graph
  qgraph(combined_coef_matrix,         
         layout = custom_layout_padded,     
         edge.labels = TRUE,             # Display edge labels
         posCol = "navy",                # Positive relationships in blue
         negCol = "red3",                # Negative relationships in red
         labels =  abbrev_names,           # Variable names
         title = plot_title,             # Plot title
         title.cex = 0.55,                  
         label.cex = 0.5,                 
         vsize = 10,
         asize = 5,
         edge.label.position = 0.55,
         edge.label.bg ="white",
         bg = bg_color,
         label.scale = FALSE,
         shape = "ellipse",
         lty = line_types,  # Use different line styles for lagged and non-lagged effects
         curve = edge_curvature)  # Apply curvature for side-by-side arrows
}


# Function to create a path diagram without lags
plot_qgraph_combined_without_lag <- function(data_ts, coef_matrix_no_lag, plot_title = "Path Diagram", bg_color = "transparent") {
  
  # Extract variable names
  abbrev_names <- colnames(data_ts)
  
  # Create an edge weight matrix initialized with zeros for combined plot
  combined_coef_matrix <- coef_matrix_no_lag
  combined_coef_matrix[,] <- 0
  
  # Store line types: default solid (1), change to dashed (2) for lagged effects
  line_types <- matrix(1, nrow = nrow(combined_coef_matrix), ncol = ncol(combined_coef_matrix))
  
  # Store curvature settings
  edge_curvature <- matrix(0, nrow = nrow(combined_coef_matrix), ncol = ncol(combined_coef_matrix))
  
  # Populate matrix with non-lagged effects (solid lines)
  nonzero_no_lag <- coef_matrix_no_lag != 0
  combined_coef_matrix[nonzero_no_lag] <- coef_matrix_no_lag[nonzero_no_lag]
  
  # Apply slight outward curvature for non-lagged edges
  edge_curvature[nonzero_no_lag] <- .85  # Slight outward curve

  # Set plot margins
  par(mar = c(0, 13, 0, 0))  
  
  # Plot the graph
  qgraph(combined_coef_matrix,         
         layout = custom_layout_padded,     
         edge.labels = TRUE,             # Display edge labels
         posCol = "navy",                # Positive relationships in blue
         negCol = "red3",                # Negative relationships in red
         labels =  abbrev_names,          # Variable names
         title = plot_title,             # Plot title
         title.cex = 0.55,                  
         label.cex = 0.5,  
         asize =5,
         edge.label.position = 0.55,
         edge.label.bg ="white",
         vsize = 10,
         bg = bg_color,
         label.scale = FALSE,
         shape = "ellipse",
         lty = line_types,  # Use different line styles for lagged and non-lagged effects
         curve = edge_curvature)  # Apply curvature for side-by-side arrows
}


#### AB Sciaenid 
coef_matrix_AB_Sciaenid_NoLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=0), "AransasBay")$coef
coef_matrix_AB_Sciaenid_YesLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=1), "AransasBay")$coef

df_Sciaenid <- data.frame(Salinity = numeric(),
                          SpottedSeatrout = numeric(),
                          PDSI = numeric(),
                          AtlanticCroaker = numeric(),
                          RedDrum = numeric(),
                          BlueCrab = numeric(),
                          stringsAsFactors = FALSE)

tiff("path_diagram_AB_Sci.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Sciaenid, coef_matrix_AB_Sciaenid_YesLag, plot_title = "Aransas Bay - Sciaenid System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Sciaenid, coef_matrix_AB_Sciaenid_NoLag, plot_title = "Aransas Bay - Sciaenid System", bg_color = "transparent")
dev.off()

#### CB Sciaenid 
coef_matrix_CB_Sciaenid_NoLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=0), "CopanoBay")$coef
coef_matrix_CB_Sciaenid_YesLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=1), "CopanoBay")$coef

tiff("path_diagram_CB_Sci.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Sciaenid, coef_matrix_CB_Sciaenid_YesLag, plot_title = "Copano Bay - Sciaenid System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Sciaenid, coef_matrix_CB_Sciaenid_NoLag, plot_title = "Copano Bay - Sciaenid System ", bg_color = "transparent")
dev.off()

#### MB Sciaenid 
coef_matrix_MB_Sciaenid_NoLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=0), "MesquiteBay")$coef
coef_matrix_MB_Sciaenid_YesLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=1), "MesquiteBay")$coef

tiff("path_diagram_MB_Sci.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Sciaenid, coef_matrix_MB_Sciaenid_YesLag, plot_title = "Mesquite Bay - Sciaenid System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Sciaenid, coef_matrix_MB_Sciaenid_NoLag, plot_title = "Mesquite Bay - Sciaenid System ", bg_color = "transparent")
dev.off()

#### GB Sciaenid 
coef_matrix_GB_Sciaenid_NoLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=0), "GalvestonBay")$coef
coef_matrix_GB_Sciaenid_YesLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=1), "GalvestonBay")$coef

tiff("path_diagram_GB_Sci.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Sciaenid, coef_matrix_GB_Sciaenid_YesLag, plot_title = "Galveston Bay - Sciaenid System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Sciaenid, coef_matrix_GB_Sciaenid_NoLag, plot_title = "Galveston Bay - Sciaenid System", bg_color = "transparent")
dev.off()

#### TB Sciaenid 
coef_matrix_TB_Sciaenid_NoLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=0), "TrinityBay")$coef
coef_matrix_TB_Sciaenid_YesLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=1), "TrinityBay")$coef

tiff("path_diagram_TB_Sci.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Sciaenid, coef_matrix_TB_Sciaenid_YesLag, plot_title = "Trinity Bay - Sciaenid System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Sciaenid, coef_matrix_TB_Sciaenid_NoLag, plot_title = "Trinity Bay - Sciaenid System", bg_color = "transparent")
dev.off()

#### EB Sciaenid 
coef_matrix_EB_Sciaenid_NoLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=0), "EastBay")$coef
coef_matrix_EB_Sciaenid_YesLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=1), "EastBay")$coef

tiff("path_diagram_EB_Sci.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Sciaenid, coef_matrix_EB_Sciaenid_YesLag, plot_title = "East Bay - Sciaenid System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Sciaenid, coef_matrix_EB_Sciaenid_NoLag, plot_title = "East Bay - Sciaenid System", bg_color = "transparent")
dev.off()

#### WB Sciaenid 
coef_matrix_WB_Sciaenid_NoLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=0), "WestBay")$coef
coef_matrix_WB_Sciaenid_YesLag <- get_part_Sciaenid_site(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=1), "WestBay")$coef

tiff("path_diagram_WB_Sci.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Sciaenid, coef_matrix_WB_Sciaenid_YesLag, plot_title = "West Bay - Sciaenid System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Sciaenid, coef_matrix_WB_Sciaenid_NoLag, plot_title = "West Bay - Sciaenid System", bg_color = "transparent")
dev.off()

#### AB keystone
coef_matrix_AB_Pred_NoLag <- get_part_Pred_site(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=0), "AransasBay")$coef
coef_matrix_AB_Pred_YesLag <- get_part_Pred_site(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=1), "AransasBay")$coef

df_Pred<- data.frame(Salinity = numeric(),
                         AlligatorGar = numeric(),
                          PDSI = numeric(),
                          Mullet = numeric(),
                          BullShark = numeric(),
                          Menhaden = numeric(),
                          stringsAsFactors = FALSE)

tiff("path_diagram_AB_Pred.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Pred, coef_matrix_AB_Pred_YesLag, plot_title = "Aransas Bay - Keystone Predator System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Pred, coef_matrix_AB_Pred_NoLag, plot_title = "Aransas Bay - Keystone Predator System", bg_color = "transparent")
dev.off()

#### CB keystone
coef_matrix_CB_Pred_NoLag <- get_part_Pred_site(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=0), "CopanoBay")$coef
coef_matrix_CB_Pred_YesLag <- get_part_Pred_site(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=1), "CopanoBay")$coef

tiff("path_diagram_CB_Pred.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Pred, coef_matrix_CB_Pred_YesLag, plot_title = "Copano Bay - Keystone Predator System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Pred, coef_matrix_CB_Pred_NoLag, plot_title = "Copano Bay - Keystone Predator System", bg_color = "transparent")
dev.off()

#### MB keystone
coef_matrix_MB_Pred_NoLag <- get_part_Pred_site(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=0), "MesquiteBay")$coef
coef_matrix_MB_Pred_YesLag <- get_part_Pred_site(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=1), "MesquiteBay")$coef

tiff("path_diagram_MB_Pred.tiff", width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Pred, coef_matrix_MB_Pred_YesLag, plot_title = "Mesquite Bay - Keystone Predator System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Pred, coef_matrix_MB_Pred_NoLag, plot_title = "Mesquite Bay - Keystone Predator System", bg_color = "transparent")
dev.off()

#### GB keystone
coef_matrix_GB_Pred_NoLag <- get_part_Pred_site_GB(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=0), "GalvestonBay")$coef
coef_matrix_GB_Pred_YesLag <- get_part_Pred_site_GB(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=1), "GalvestonBay")$coef

tiff("path_diagram_GB_Pred.tiff",width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Pred, coef_matrix_GB_Pred_YesLag, plot_title = "Galveston Bay - Keystone Predator System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Pred, coef_matrix_GB_Pred_NoLag, plot_title = "Galveston Bay - Keystone Predator System", bg_color = "transparent")
dev.off()

#### TB keystone
coef_matrix_TB_Pred_NoLag <- get_part_Pred_site_GB(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=0), "TrinityBay")$coef
coef_matrix_TB_Pred_YesLag <- get_part_Pred_site_GB(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=1), "TrinityBay")$coef

tiff("path_diagram_TB_Pred.tiff",width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Pred, coef_matrix_TB_Pred_YesLag, plot_title = "Trinity Bay - Keystone Predator System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Pred, coef_matrix_TB_Pred_NoLag, plot_title = "Trinity Bay - Keystone Predator System", bg_color = "transparent")
dev.off()

#### EB keystone
coef_matrix_EB_Pred_NoLag <- get_part_Pred_site_GB(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=0), "EastBay")$coef
coef_matrix_EB_Pred_YesLag <- get_part_Pred_site_GB(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=1), "EastBay")$coef

tiff("path_diagram_EB_Pred.tiff",width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Pred, coef_matrix_EB_Pred_YesLag, plot_title = "East Bay - Keystone Predator System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Pred, coef_matrix_EB_Pred_NoLag, plot_title = "East Bay - Keystone Predator System", bg_color = "transparent")
dev.off()

#### WB keystone
coef_matrix_WB_Pred_NoLag <- get_part_Pred_site_GB(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=0), "WestBay")$coef
coef_matrix_WB_Pred_YesLag <- get_part_Pred_site_GB(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=1), "WestBay")$coef

tiff("path_diagram_WB_Pred.tiff",width = 55,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Pred, coef_matrix_WB_Pred_YesLag, plot_title = "West Bay - Keystone Predator System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Pred, coef_matrix_WB_Pred_NoLag, plot_title = "West Bay - Keystone Predator System", bg_color = "transparent")
dev.off()


###### make 'empty' path diagrams as conceptual models 


# duplicate  original sciaenid matrix and  fill in with dummy values
coef_matrix_dummy_Sciaenid <- matrix(1, 
                                     nrow = nrow(coef_matrix_AB_Sciaenid_NoLag), 
                                     ncol = ncol(coef_matrix_AB_Sciaenid_NoLag), 
                                     dimnames = dimnames(coef_matrix_AB_Sciaenid_NoLag))
coef_matrix_dummy_Pred <- matrix(1, 
                                 nrow = nrow(coef_matrix_AB_Pred_NoLag), 
                                 ncol = ncol(coef_matrix_AB_Pred_NoLag), 
                                 dimnames = dimnames(coef_matrix_AB_Pred_NoLag))


# define custom layout with (x, y) coordinates for path diagrams
custom_layout <- matrix(c(
  -1,  1.9, # salinity
  1,  1,  # spottedseatrout
  -1,  1.1,   # PDSI
  0, 1.7,   # croaker
  1,  2,   # drum 
  0, 1.3    # blue crab
), ncol = 2, byrow = TRUE)

# function to make conceptual path diagram
plot_qgraph_panel_conceptual <- function(data_ts, coef_matrix, plot_title = "Fitted DSEM Path Diagram") {
  abbrev_names <- colnames(data_ts)
  if (is.null(abbrev_names)) {
    stop("The time series object does not have column names.")
  }
  par(mar = c(25, 25, 25, 25))  
  num_vars <- ncol(coef_matrix)
  edge_curvature <- matrix(0.8, nrow = num_vars, ncol = num_vars) 
  edge_colors <- matrix("black", nrow = num_vars, ncol = num_vars) 
  qgraph(coef_matrix,         
         layout = custom_layout,     
         edge.labels = FALSE,              
         posCol = "black",                 
         negCol = "black",                 
         labels = abbrev_names,            
         title = plot_title,              
         title.cex = 0.7,                  
         label.cex = 0.5,                 
         vsize = 10,
         asize = 3,
         edge.label.position = 0.55,
         edge.label.bg = "white",
         label.scale = FALSE,
         shape = "ellipse",                
         directed = TRUE,  
         arrows = TRUE,   
         edge.color = edge_colors, 
         edge.width = 2.5,  
         curve = edge_curvature)
}

# replacing some 1s with 0s to 'remove' non modeled relationships 
coef_matrix_dummy_Sciaenid["Salinity", "PDSI"] <- 0
coef_matrix_dummy_Sciaenid["SpottedSeatrout", "PDSI"] <- 0
coef_matrix_dummy_Sciaenid["RedDrum", "PDSI"] <- 0
coef_matrix_dummy_Sciaenid["BlueCrabSmall", "PDSI"] <- 0
coef_matrix_dummy_Sciaenid["Atlanticcroaker", "PDSI"] <- 0
coef_matrix_dummy_Sciaenid["SpottedSeatrout", "Salinity"] <- 0
coef_matrix_dummy_Sciaenid["RedDrum", "Salinity"] <- 0
coef_matrix_dummy_Sciaenid["BlueCrabSmall", "Salinity"] <- 0
coef_matrix_dummy_Sciaenid["Atlanticcroaker", "Salinity"] <- 0
coef_matrix_dummy_Sciaenid["SpottedSeatrout", "RedDrum"] <- 0
coef_matrix_dummy_Sciaenid["RedDrum", "SpottedSeatrout"] <- 0
coef_matrix_dummy_Sciaenid["BlueCrabSmall", "Atlanticcroaker"] <- 0
coef_matrix_dummy_Sciaenid["Atlanticcroaker", "BlueCrabSmall"] <- 0

coef_matrix_dummy_Pred["Salinity", "PDSI"] <- 0
coef_matrix_dummy_Pred["BullShark", "PDSI"] <- 0
coef_matrix_dummy_Pred["AlligatorGar", "PDSI"] <- 0
coef_matrix_dummy_Pred["Mullet", "PDSI"] <- 0
coef_matrix_dummy_Pred["Menhaden", "PDSI"] <- 0
coef_matrix_dummy_Pred["BullShark", "Salinity"] <- 0
coef_matrix_dummy_Pred["AlligatorGar", "Salinity"] <- 0
coef_matrix_dummy_Pred["Mullet", "Salinity"] <- 0
coef_matrix_dummy_Pred["Menhaden", "Salinity"] <- 0
coef_matrix_dummy_Pred["BullShark", "AlligatorGar"] <- 0
coef_matrix_dummy_Pred["AlligatorGar", "BullShark"] <- 0
coef_matrix_dummy_Pred["Mullet", "Menhaden"] <- 0
coef_matrix_dummy_Pred["Menhaden", "Mullet"] <- 0


# make and save plots
tiff("~/Desktop/SciaenidSystem.tiff",width = 85,  height = 110, units = "mm", res = 300,bg = "transparent")
par(mfrow = c(1, 1))
plot_qgraph_panel_conceptual(df_Sciaenid, coef_matrix_dummy_Sciaenid, plot_title = "Sciaenid System")
dev.off()

tiff("~/Desktop/KeystoneSystem.tiff",width = 85,  height = 110, units = "mm", res = 300,bg = "transparent")
par(mfrow = c(1, 1))
plot_qgraph_panel_conceptual(df_Pred, coef_matrix_dummy_Pred, plot_title = "Keystone Predator System")
dev.off()


###### make synthesis plot that shows all significant paths from all lowest AIC models


# first extract data/info from signficant paths
extract_sig_paths <- function(model, model_name) {
  df <- summary(model) %>%
    as.data.frame() %>%
    mutate( model = model_name,
      
      # extract bay 
      minor_bay = if_else(
        str_detect(second, "_"),
        str_extract(second, "(?<=_).+$"),
        NA_character_),
      
      # extract response variable 
      response = if_else(
        str_detect(second, "_"),
        str_extract(second, "^[^_]+"),
        second),
      
      # extract predictor variable 
      predictor = if_else(
        str_detect(first, "_"),
        str_extract(first, "^[^_]+"),
        first),
      
      # set Estimate to 0 if not significant
      Estimate = if_else(p_value >= 0.05, 0, Estimate)) %>%
    
    # remove lag0 self-relationships
    filter(!(predictor == response & lag == 0)) %>%
    
    # keep only essentials
    select(predictor, response, minor_bay,lag,Estimate,p_value,model)
  return(df)
}

# apply function to lowest AIC models and bind data
sig_AB_Sciaenid <- extract_sig_paths(fit_semAB_Sciaenid_notrophics, "AB_Sciaenid_notrophics")
sig_GB_Sciaenid  <- extract_sig_paths(fit_semGB_Sciaenid_fulltopdown, "GB_Sciaenid_fulltopdown")
sig_AB_Pred      <- extract_sig_paths(fit_semAB_Pred_fulltopdown,    "AB_Pred_fulltopdown")
sig_GB_Pred      <- extract_sig_paths(fit_semGB_Pred_fullbottomup,   "GB_Pred_fullbottomup")

sig_all <- bind_rows(sig_AB_Sciaenid, sig_GB_Sciaenid,sig_AB_Pred, sig_GB_Pred)

# a whole bunch of inefficient reformatting code to make plotting smoother/prettier
sig_all2 <- sig_all %>%
  mutate(context = paste0(minor_bay, "-", lag),
    sign = case_when(
      Estimate > 0 ~ "positive",
      Estimate < 0 ~ "negative",
      TRUE ~ "none"),
    pair_label = paste0(predictor, " → ", response))

format_context_label <- function(ctx) {
  ctx <- as.character(ctx)
  na_idx <- is.na(ctx)
  out <- rep(NA_character_, length(ctx))
  m <- str_match(ctx, "^([A-Za-z]+Bay)-([0-9]+)$")
  bay_raw <- m[,2]
  lag <- m[,3]
    bay_spaced <- ifelse(
    is.na(bay_raw),
    NA_character_,
    str_replace_all(bay_raw, "(?<=.)([A-Z])", " \\1"))
  bay_spaced <- str_replace(bay_spaced, "Bay$", "Bay") 
  out[!is.na(bay_raw)] <- paste0(bay_spaced[!is.na(bay_raw)], " - ", lag[!is.na(bay_raw)], " yr lag")
  out[na_idx] <- NA_character_
  return(out)
}

sig_all2 <- sig_all2 %>%
  mutate(context = paste0(minor_bay, "-", lag))

sig_all2 <- sig_all2 %>%
  mutate(pair_label = str_replace_all(pair_label, c("AllMullet" = "Mullet", 
                                               "BlueCrabSmall" = "BlueCrab", 
                                               "Atlanticcroaker" = "AtlanticCroaker", 
                                               "AllMenhaden" = "Menhaden")))

bay_order <- c("AransasBay", "MesquiteBay", "CopanoBay", "WestBay",
               "GalvestonBay", "EastBay", "TrinityBay")

lag_order <- c("0", "1")

context_order <- paste0(rep(bay_order, each = 2), "-", lag_order)

context_label_levels <- format_context_label(context_order)

sig_all2 <- sig_all2 %>%
  mutate(context_label = format_context_label(context))

all_pairs <- unique(sig_all2$pair_label)

# customizing order of y axis (path) labels 
order_PDSI <- grep("^PDSI", all_pairs, value = TRUE)

order_Salinity <- grep("^Salinity", all_pairs, value = TRUE)

auto_vars <- c("Mullet", "Menhaden", "RedDrum",
               "SpottedSeatrout", "BullShark", "AlligatorGar")

order_AR <- unlist(lapply(auto_vars, function(v) {
  grep(paste0("^", v, " → ", v, "$"), all_pairs, value = TRUE)
}))

order_MM <- grep("^(Mullet|Menhaden)", all_pairs, value = TRUE)
order_MM <- setdiff(order_MM, order_AR)

top_pred_vars <- c("RedDrum", "SpottedSeatrout", "AlligatorGar", "BullShark")

order_top_preds <- unlist(lapply(top_pred_vars, function(v) {
  setdiff(grep(paste0("^", v), all_pairs, value = TRUE), order_AR)
}))

new_pair_order <- c(order_PDSI,order_Salinity,order_AR,order_MM,order_top_preds)

sig_all2 <- sig_all2 %>% 
  mutate(pair_label = factor(pair_label, levels = new_pair_order))


# build full_grid with all formatting and plotting specs
full_grid <- expand_grid(
  pair_label = unique(sig_all2$pair_label),
  context = context_order) %>%
  left_join(sig_all2 %>% select(pair_label, context, sign, context_label),
    by = c("pair_label", "context")) %>%
  mutate(
    has_relationship = !is.na(sign),
    sign = ifelse(is.na(sign), "none", sign),
    sign = factor(sign, levels = c("positive", "negative", "none")),
    pair_label = factor(pair_label, levels = new_pair_order))

full_grid <- full_grid %>%
  mutate(context_label = factor(format_context_label(context), levels = context_label_levels))

# specifying placement of lines to make plot prettier
vline_pos <- which(context_label_levels == format_context_label("CopanoBay-1")) + 0.5

y_levels <- levels(full_grid$pair_label)

y_pos <- function(label) {
  which(y_levels == label)}

lines_to_add <- c(
  y_pos("PDSI → AlligatorGar") + 0.5,
  y_pos("Salinity → AlligatorGar") + 0.5,
  y_pos("AlligatorGar → AlligatorGar") + 0.5,
  y_pos("Menhaden → AlligatorGar") + 0.5)


# make dot plot / heat plot to show sig path coefficients for diff paths and diff bays
dotplot <- ggplot(full_grid, aes(x = context_label, y = pair_label)) +
  geom_tile(fill = "white", color = "black") +
  geom_vline(xintercept = vline_pos, color = "black",linewidth = 0.5) +
  geom_hline(yintercept = lines_to_add,linetype = "dashed",
    color = "grey40",linewidth = 0.5)+
  geom_point(data = subset(full_grid, has_relationship),
    aes(color = sign, fill = sign),size = 3.5,shape = 21,stroke = 0.1,alpha = 1.0) +
  scale_color_manual(
    values = c("positive" = "black",
               "negative" = "black",
               "none" = "black"),drop = FALSE) +
  scale_fill_manual(
    values = c("positive" = "dodgerblue1",
               "negative" = "firebrick2",
               "none" = "white"),drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "none") +
  labs(x = "Bay (and lag)", y = "Path")+
  annotate("text", x = vline_pos / 2, y = length(levels(full_grid$pair_label)),
    label = "Mission-Aransas Estuary",fontface = "bold",size = 3.5,vjust = -2) +
  annotate("text",x = (vline_pos + length(context_label_levels)) / 2, y = length(levels(full_grid$pair_label)),
    label = "Trinity-San Jacinto Estuary",fontface = "bold",size = 3.5,vjust = -2) +
  geom_text( data = data.frame(
      label = c("PDSI Effects", "Salinity Effects", "Top-Down Effects"),
      y = c(7.5, 16, 34.5), x = length(context_label_levels) + 0.5),
    aes(x = x, y = y, label = label), hjust = 0, vjust = -1.5,angle = 270, size = 2.5) +
  geom_text(data = data.frame(label = "Density Dependent", y = 23.5,x = length(context_label_levels) + 0.5),
    aes(x = x, y = y, label = label), hjust = 0, vjust = -2.5, angle = 270, size = 2.5) +
  geom_text(data = data.frame(label = "Effects",y = 21.5, x = length(context_label_levels) + 0.55),
    aes(x = x, y = y, label = label), hjust = 0, vjust = -0.7,angle = 270, size = 2.5)+
  geom_text(data = data.frame(label = "Bottom-Up", y = 27.5,x = length(context_label_levels) + 0.5),
            aes(x = x, y = y, label = label), hjust = 0, vjust = -2.5, angle = 270, size = 2.5) +
  geom_text(data = data.frame(label = "Effects",y = 26.5, x = length(context_label_levels) + 0.55),
            aes(x = x, y = y, label = label), hjust = 0, vjust = -0.7,angle = 270, size = 2.5)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 20, r = 30, b = 10, l = 10))
dotplot


ggsave("dsemdotplot.tiff", dotplot, width = 170, height = 170, dpi = 300, units = "mm", bg = "white")

