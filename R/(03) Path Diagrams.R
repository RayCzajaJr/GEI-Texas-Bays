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


# Extract coefficients
get_part_Sciaenid = function(x){
  vars = c("Salinity","SpottedSeatrout","PDSI","Atlanticcroaker","RedDrum", "BlueCrabSmall")
  index = sapply( vars, FUN=\(y) grep(y,rownames(x$coef))[1] )
  x$coef = x$coef[index,index]
  dimnames(x$coef) = list( vars, vars )
  return(x)
}

get_part_Pred = function(x){
  vars = c("Salinity","AlligatorGar","PDSI","Mullet","BullShark", "Menhaden")
  index = sapply( vars, FUN=\(y) grep(y,rownames(x$coef))[1] )
  x$coef = x$coef[index,index]
  dimnames(x$coef) = list( vars, vars )
  return(x)
}


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
  par(mar = c(0, 0, 0, 0))  
  
  # Plot the graph
  qgraph(combined_coef_matrix,         
         layout = custom_layout_padded,     
         edge.labels = TRUE,             # Display edge labels
         posCol = "navy",                # Positive relationships in blue
         negCol = "red3",                # Negative relationships in red
         labels =  abbrev_names,           # Variable names
         title = plot_title,             # Plot title
         title.cex = 0.7,                  
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
  edge_curvature[nonzero_no_lag] <- 1.5  # Slight outward curve

  
  # Set plot margins
  par(mar = c(0, 0, 0, 0))  
  
  # Plot the graph
  qgraph(combined_coef_matrix,         
         layout = custom_layout_padded,     
         edge.labels = TRUE,             # Display edge labels
         posCol = "navy",                # Positive relationships in blue
         negCol = "red3",                # Negative relationships in red
         labels =  abbrev_names,          # Variable names
         title = plot_title,             # Plot title
         title.cex = 0.7,                  
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


coef_matrix_AB_Sciaenid_NoLag <- get_part_Sciaenid(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=0))$coef
coef_matrix_AB_Sciaenid_YesLag <- get_part_Sciaenid(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=1))$coef

df_Sciaenid <- data.frame(Salinity = numeric(),
                          SpottedSeatrout = numeric(),
                          PDSI = numeric(),
                          AtlanticCroaker = numeric(),
                          RedDrum = numeric(),
                          BlueCrab = numeric(),
                          stringsAsFactors = FALSE)

tiff("path_diagram_AB_Sci.tiff", width = 85,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Sciaenid, coef_matrix_AB_Sciaenid_YesLag, plot_title = "Aransas Bay - Sciaenid System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Sciaenid, coef_matrix_AB_Sciaenid_NoLag, plot_title = "Aransas Bay - Sciaenid System", bg_color = "transparent")
dev.off()


#### GB Sciaenid 


coef_matrix_GB_Sciaenid_NoLag <- get_part_Sciaenid(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=0))$coef
coef_matrix_GB_Sciaenid_YesLag <- get_part_Sciaenid(as_fitted_DAG(fit_semGB_Sciaenid_fulltopdown, lag=1))$coef

df_Sciaenid <- data.frame(Salinity = numeric(),
                          SpottedSeatrout = numeric(),
                          PDSI = numeric(),
                          AtlanticCroaker = numeric(),
                          RedDrum = numeric(),
                          BlueCrab = numeric(),
                          stringsAsFactors = FALSE)

tiff("path_diagram_GB_Sci.tiff", width = 85,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Sciaenid, coef_matrix_GB_Sciaenid_YesLag, plot_title = "Galveston Bay - Sciaenid System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Sciaenid, coef_matrix_GB_Sciaenid_NoLag, plot_title = "Galveston Bay - Sciaenid System", bg_color = "transparent")
dev.off()

#### AB keystone


coef_matrix_AB_Pred_NoLag <- get_part_Pred(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=0))$coef
coef_matrix_AB_Pred_YesLag <- get_part_Pred(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=1))$coef

df_Pred<- data.frame(Salinity = numeric(),
                         AlligatorGar = numeric(),
                          PDSI = numeric(),
                          Mullet = numeric(),
                          BullShark = numeric(),
                          Menhaden = numeric(),
                          stringsAsFactors = FALSE)

tiff("path_diagram_AB_Pred.tiff", width = 85,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Pred, coef_matrix_AB_Pred_YesLag, plot_title = "Aransas Bay - Keystone Predator System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Pred, coef_matrix_AB_Pred_NoLag, plot_title = "Aransas Bay - Keystone Predator System", bg_color = "transparent")
dev.off()


#### GB keystone


coef_matrix_GB_Pred_NoLag <- get_part_Pred(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=0))$coef
coef_matrix_GB_Pred_YesLag <- get_part_Pred(as_fitted_DAG(fit_semGB_Pred_fullbottomup, lag=1))$coef

df_Pred<- data.frame(Salinity = numeric(),
                     AlligatorGar = numeric(),
                     PDSI = numeric(),
                     Mullet = numeric(),
                     BullShark = numeric(),
                     Menhaden = numeric(),
                     stringsAsFactors = FALSE)

tiff("path_diagram_GB_Pred.tiff",width = 85,  height = 110, units = "mm", res = 300,bg = "transparent")
plot_qgraph_combined_with_lag(df_Pred, coef_matrix_GB_Pred_YesLag, plot_title = "Galveston Bay - Keystone Predator System", bg_color = "white")
plot_qgraph_combined_without_lag(df_Pred, coef_matrix_GB_Pred_NoLag, plot_title = "Galveston Bay - Keystone Predator System", bg_color = "transparent")
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
