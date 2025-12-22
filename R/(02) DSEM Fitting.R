
library(dsem)
library(ggplot2)
library(ggpubr)
library(ggraph)
library(phylopath)
library(dplyr)
library(ggdag)
library(readxl)
library(writexl)
library(ggraph)
library(qgraph)
library(DHARMa)


########## Make six DSEMS per trophic model and major bay, based on the 6 different scenarios 
########## Determine the lowest AIC model per trophic model and major bay


# First System and Major Bay: Sciaenid Trophic System and Aransas Bay 

# get data into a time series object from all numeric columns except YEAR
AB_Sciaenid_TS <- ts(
  AransasBay_Sciaenid_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(AransasBay_Sciaenid_Wide_Trans$YEAR)),  
  end = c(max(AransasBay_Sciaenid_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

# specify the dsem model for the sciaenid system for aransas bay for zero and 1 year lags

# full bottom up
semAB_Sciaenid_fullbottomup <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Prey predicting Sciaenid
Atlanticcroaker_MesquiteBay -> RedDrum_MesquiteBay, 0, b1
Atlanticcroaker_AransasBay -> RedDrum_AransasBay, 0, b2
Atlanticcroaker_CopanoBay -> RedDrum_CopanoBay, 0, b3

BlueCrabSmall_MesquiteBay -> RedDrum_MesquiteBay, 0, b4
BlueCrabSmall_AransasBay -> RedDrum_AransasBay, 0, b5
BlueCrabSmall_CopanoBay -> RedDrum_CopanoBay, 0, b6

Atlanticcroaker_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, b7
Atlanticcroaker_AransasBay -> SpottedSeatrout_AransasBay, 0, b8
Atlanticcroaker_CopanoBay -> SpottedSeatrout_CopanoBay, 0, b9

BlueCrabSmall_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, b10
BlueCrabSmall_AransasBay -> SpottedSeatrout_AransasBay, 0, b11
BlueCrabSmall_CopanoBay -> SpottedSeatrout_CopanoBay, 0, b12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12

RedDrum_MesquiteBay -> RedDrum_MesquiteBay, 1, l25
RedDrum_AransasBay -> RedDrum_AransasBay, 1, l26
RedDrum_CopanoBay -> RedDrum_CopanoBay, 1, l27

SpottedSeatrout_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l28
SpottedSeatrout_AransasBay -> SpottedSeatrout_AransasBay, 1, l29
SpottedSeatrout_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l30
"

fit_semAB_Sciaenid_fullbottomup = dsem(sem = semAB_Sciaenid_fullbottomup,
                          tsdata = AB_Sciaenid_TS,
                          control = dsem_control(
                            quiet = TRUE))
summary(fit_semAB_Sciaenid_fullbottomup)


#### Preliminary inspection of some path diagrams.... the next script will more fully develop diagrams for all trophic models and bays 

# extract the coefficients
get_part_Sciaenid = function(x){
  vars = c("Salinity","SpottedSeatrout","PDSI","Atlanticcroaker","RedDrum", "BlueCrabSmall")
  index = sapply( vars, FUN=\(y) grep(y,rownames(x$coef))[1] )
  x$coef = x$coef[index,index]
  dimnames(x$coef) = list( vars, vars )
  return(x)
}

# function to make path diagram
plot_qgraph_panel <- function(data_ts, coef_matrix, plot_title = "Fitted DSEM Path Diagram") {
  
  # Extract the column names (variable names) from the time series object
  abbrev_names <- colnames(data_ts)
  
  # Check if the time series object has column names
  if (is.null(abbrev_names)) {
    stop("The time series object does not have column names.")
  }
  
  # Set plot margins (bottom, left, top, right) to add cushion space
  par(mar = c(25, 25, 25, 25))  
  
  # Plot the graph
  qgraph(coef_matrix,         
         layout = "groups",              # Layout style
         edge.labels = TRUE,             # Display edge labels
         posCol = "navy",                # Set positive relationships to blue
         negCol = "red3",                 # Optionally set negative relationships to red
         labels = abbrev_names,          # Use the automatically extracted variable names
         title = plot_title,             # Use the specified plot title
         title.cex = 1,                  # Increase font size of the title
         label.cex = .8,                  # Increase font size of the variable names
         vsize = 10,
         label.scale = FALSE,
         shape = "ellipse")
}


# get the coefficient into a matrix for the path diagrams
coef_matrix_AB_Sciaenid_NoLag <- get_part_Sciaenid(as_fitted_DAG(fit_semAB_Sciaenid_fullbottomup, lag=0))
coef_matrix_AB_Sciaenid_YesLag <- get_part_Sciaenid(as_fitted_DAG(fit_semAB_Sciaenid_fullbottomup, lag=1))

# create an empty data frame with specified column names
df_Sciaenid <- data.frame(Salinity = numeric(),
                          SpottedSeatrout = numeric(),
                          PDSI = numeric(),
                          Atlanticcroaker = numeric(),
                          RedDrum = numeric(),
                          BlueCrabSmall = numeric(),
                          stringsAsFactors = FALSE)

png("~/Desktop/AB_Sciaenid.png", width = 1200, height = 2000, res = 200)
par(mfrow = c(2, 1))
plot_qgraph_panel(df_Sciaenid, coef_matrix_AB_Sciaenid_NoLag$coef,  plot_title = "Sciaenid-AB-No Lag")
plot_qgraph_panel(df_Sciaenid, coef_matrix_AB_Sciaenid_YesLag$coef, plot_title = "Sciaenid-AB-1 Year Lag")
par(mfrow = c(1, 1))
dev.off()

# no DD bottom up model
semAB_Sciaenid_bottomup_noDD <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Prey predicting Sciaenid
Atlanticcroaker_MesquiteBay -> RedDrum_MesquiteBay, 0, b1
Atlanticcroaker_AransasBay -> RedDrum_AransasBay, 0, b2
Atlanticcroaker_CopanoBay -> RedDrum_CopanoBay, 0, b3

BlueCrabSmall_MesquiteBay -> RedDrum_MesquiteBay, 0, b4
BlueCrabSmall_AransasBay -> RedDrum_AransasBay, 0, b5
BlueCrabSmall_CopanoBay -> RedDrum_CopanoBay, 0, b6

Atlanticcroaker_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, b7
Atlanticcroaker_AransasBay -> SpottedSeatrout_AransasBay, 0, b8
Atlanticcroaker_CopanoBay -> SpottedSeatrout_CopanoBay, 0, b9

BlueCrabSmall_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, b10
BlueCrabSmall_AransasBay -> SpottedSeatrout_AransasBay, 0, b11
BlueCrabSmall_CopanoBay -> SpottedSeatrout_CopanoBay, 0, b12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12

"

fit_semAB_Sciaenid_bottomup_noDD = dsem(sem = semAB_Sciaenid_bottomup_noDD,
                          tsdata = AB_Sciaenid_TS,
                          control = dsem_control(
                            quiet = TRUE))

summary(fit_semAB_Sciaenid_bottomup_noDD)

# abiotic only model
semAB_Sciaenid_abiotic <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12


# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12
"

fit_semAB_Sciaenid_abiotic = dsem(sem = semAB_Sciaenid_abiotic,
                          tsdata = AB_Sciaenid_TS,
                          control = dsem_control(
                            quiet = TRUE))

summary(fit_semAB_Sciaenid_abiotic)

# Full top down model
semAB_Sciaenid_fulltopdown <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Sciaenid predicting prey
RedDrum_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, b1
RedDrum_AransasBay -> Atlanticcroaker_AransasBay, 0, b2
RedDrum_CopanoBay -> Atlanticcroaker_CopanoBay, 0, b3

RedDrum_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, b4
RedDrum_AransasBay -> BlueCrabSmall_AransasBay, 0, b5
RedDrum_CopanoBay -> BlueCrabSmall_CopanoBay, 0, b6

SpottedSeatrout_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, b7
SpottedSeatrout_AransasBay -> Atlanticcroaker_AransasBay, 0, b8
SpottedSeatrout_CopanoBay -> Atlanticcroaker_CopanoBay, 0, b9

SpottedSeatrout_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, b10
SpottedSeatrout_AransasBay -> BlueCrabSmall_AransasBay, 0, b11
SpottedSeatrout_CopanoBay -> BlueCrabSmall_CopanoBay, 0, b12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12

RedDrum_MesquiteBay -> RedDrum_MesquiteBay, 1, l25
RedDrum_AransasBay -> RedDrum_AransasBay, 1, l26
RedDrum_CopanoBay -> RedDrum_CopanoBay, 1, l27

SpottedSeatrout_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l28
SpottedSeatrout_AransasBay -> SpottedSeatrout_AransasBay, 1, l29
SpottedSeatrout_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l30
"

# Fit the model for prey prediction in Aransas Bay
fit_semAB_Sciaenid_fulltopdown = dsem(sem = semAB_Sciaenid_fulltopdown,
                          tsdata = AB_Sciaenid_TS,
                          control = dsem_control(
                            quiet = TRUE))

summary(fit_semAB_Sciaenid_fulltopdown)



# top down no DD model
semAB_Sciaenid_topdown_noDD <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Sciaenid predicting prey
RedDrum_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, b1
RedDrum_AransasBay -> Atlanticcroaker_AransasBay, 0, b2
RedDrum_CopanoBay -> Atlanticcroaker_CopanoBay, 0, b3

RedDrum_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, b4
RedDrum_AransasBay -> BlueCrabSmall_AransasBay, 0, b5
RedDrum_CopanoBay -> BlueCrabSmall_CopanoBay, 0, b6

SpottedSeatrout_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, b7
SpottedSeatrout_AransasBay -> Atlanticcroaker_AransasBay, 0, b8
SpottedSeatrout_CopanoBay -> Atlanticcroaker_CopanoBay, 0, b9

SpottedSeatrout_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, b10
SpottedSeatrout_AransasBay -> BlueCrabSmall_AransasBay, 0, b11
SpottedSeatrout_CopanoBay -> BlueCrabSmall_CopanoBay, 0, b12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12
"

# Fit the model for prey prediction in Aransas Bay
fit_semAB_Sciaenid_topdown_noDD  = dsem(sem = semAB_Sciaenid_topdown_noDD,
                               tsdata = AB_Sciaenid_TS,
                               control = dsem_control(
                                 quiet = TRUE))
summary(fit_semAB_Sciaenid_topdown_noDD)


# No Prey Model
semAB_Sciaenid_notrophics <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12

RedDrum_MesquiteBay -> RedDrum_MesquiteBay, 1, l25
RedDrum_AransasBay -> RedDrum_AransasBay, 1, l26
RedDrum_CopanoBay -> RedDrum_CopanoBay, 1, l27

SpottedSeatrout_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l28
SpottedSeatrout_AransasBay -> SpottedSeatrout_AransasBay, 1, l29
SpottedSeatrout_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l30
"

# Fit the No Prey  model
fit_semAB_Sciaenid_notrophics = dsem(sem = semAB_Sciaenid_notrophics,
                              tsdata = AB_Sciaenid_TS,
                              control = dsem_control(
                                quiet = TRUE))

summary(fit_semAB_Sciaenid_notrophics)

# Compare AICs of all six models
AIC(fit_semAB_Sciaenid_abiotic)
AIC(fit_semAB_Sciaenid_notrophics) #WINNER
AIC(fit_semAB_Sciaenid_bottomup_noDD)
AIC(fit_semAB_Sciaenid_topdown_noDD)
AIC(fit_semAB_Sciaenid_fullbottomup)
AIC(fit_semAB_Sciaenid_fulltopdown)

# Export model parameters (eg p values) as excel file
table_semAB_Sciaenid<-summary(fit_semAB_Sciaenid_notrophics)
write_xlsx(table_semAB_Sciaenid, path = "table_semAB_Sciaenid.xlsx")



###################################################################




# Repeating the above but for the...

# Second System and Major Bay: Keystone Predator Trophic System and Aransas Bay 

#get data into a time series object from all numeric columns except YEAR
AB_Pred_TS <- ts(
  AransasBay_Pred_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(AransasBay_Pred_Wide_Trans$YEAR)),  
  end = c(max(AransasBay_Pred_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

# Full bottom up model
semAB_Pred_fullbottomup <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Prey predicting predators
AllMullet_MesquiteBay -> BullShark_MesquiteBay, 0, b1
AllMullet_AransasBay -> BullShark_AransasBay, 0, b2
AllMullet_CopanoBay -> BullShark_CopanoBay, 0, b3

AllMenhaden_MesquiteBay -> BullShark_MesquiteBay, 0, b4
AllMenhaden_AransasBay -> BullShark_AransasBay, 0, b5
AllMenhaden_CopanoBay -> BullShark_CopanoBay, 0, b6

AllMullet_MesquiteBay -> AlligatorGar_MesquiteBay, 0, b7
AllMullet_AransasBay -> AlligatorGar_AransasBay, 0, b8
AllMullet_CopanoBay -> AlligatorGar_CopanoBay, 0, b9

AllMenhaden_MesquiteBay -> AlligatorGar_MesquiteBay, 0, b10
AllMenhaden_AransasBay -> AlligatorGar_AransasBay, 0, b11
AllMenhaden_CopanoBay -> AlligatorGar_CopanoBay, 0, b12

# Lagged relationships (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

AllMullet_MesquiteBay -> AllMullet_MesquiteBay, 1, l13
AllMullet_AransasBay -> AllMullet_AransasBay, 1, l14
AllMullet_CopanoBay -> AllMullet_CopanoBay, 1, l15

AllMenhaden_MesquiteBay -> AllMenhaden_MesquiteBay, 1, l16
AllMenhaden_AransasBay -> AllMenhaden_AransasBay, 1, l17
AllMenhaden_CopanoBay -> AllMenhaden_CopanoBay, 1, l18

BullShark_MesquiteBay -> BullShark_MesquiteBay, 1, l25
BullShark_AransasBay -> BullShark_AransasBay, 1, l26
BullShark_CopanoBay -> BullShark_CopanoBay, 1, l27

AlligatorGar_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l28
AlligatorGar_AransasBay -> AlligatorGar_AransasBay, 1, l29
AlligatorGar_CopanoBay -> AlligatorGar_CopanoBay, 1, l30
"

fit_semAB_Pred_fullbottomup = dsem(sem = semAB_Pred_fullbottomup,
                      tsdata = AB_Pred_TS,
                      control = dsem_control(
                        quiet = TRUE))
summary(fit_semAB_Pred_fullbottomup)


# Abiotic Only Model
semAB_Pred_abiotic <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Lagged relationships for predators (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12
"

fit_semAB_Pred_abiotic = dsem(sem = semAB_Pred_abiotic,
                                  tsdata = AB_Pred_TS,
                                  control = dsem_control(
                                    quiet = TRUE))
summary(fit_semAB_Pred_abiotic)


# No Density Dependence bottom up Model
semAB_Pred_bottomupnoDD <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Prey predicting predators
AllMullet_MesquiteBay -> BullShark_MesquiteBay, 0, b1
AllMullet_AransasBay -> BullShark_AransasBay, 0, b2
AllMullet_CopanoBay -> BullShark_CopanoBay, 0, b3

AllMenhaden_MesquiteBay -> BullShark_MesquiteBay, 0, b4
AllMenhaden_AransasBay -> BullShark_AransasBay, 0, b5
AllMenhaden_CopanoBay -> BullShark_CopanoBay, 0, b6

AllMullet_MesquiteBay -> AlligatorGar_MesquiteBay, 0, b7
AllMullet_AransasBay -> AlligatorGar_AransasBay, 0, b8
AllMullet_CopanoBay -> AlligatorGar_CopanoBay, 0, b9

AllMenhaden_MesquiteBay -> AlligatorGar_MesquiteBay, 0, b10
AllMenhaden_AransasBay -> AlligatorGar_AransasBay, 0, b11
AllMenhaden_CopanoBay -> AlligatorGar_CopanoBay, 0, b12

# Lagged relationships for predators (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888
"

# Fit the No Density Dependence bottom up  model
fit_semAB_Pred_bottomupnoDD = dsem(sem = semAB_Pred_bottomupnoDD,
                        tsdata = AB_Pred_TS,
                        control = dsem_control(
                          quiet = TRUE))
summary(fit_semAB_Pred_bottomupnoDD)


# No trophics Model
semAB_Pred_notrophics <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Lagged relationships for predators (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

BullShark_MesquiteBay -> BullShark_MesquiteBay, 1, l25
BullShark_AransasBay -> BullShark_AransasBay, 1, l26
BullShark_CopanoBay -> BullShark_CopanoBay, 1, l27

AlligatorGar_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l28
AlligatorGar_AransasBay -> AlligatorGar_AransasBay, 1, l29
AlligatorGar_CopanoBay -> AlligatorGar_CopanoBay, 1, l30
"

# Fit the No trophics SEM model
fit_semAB_Pred_notrophics= dsem(sem = semAB_Pred_notrophics,
                          tsdata = AB_Pred_TS,
                          control = dsem_control(
                            quiet = TRUE))
summary(fit_semAB_Pred_notrophics)


# Full model for top down
semAB_Pred_fulltopdown <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Predators predicting prey
BullShark_MesquiteBay -> AllMullet_MesquiteBay, 0, b1
BullShark_AransasBay -> AllMullet_AransasBay, 0, b2
BullShark_CopanoBay -> AllMullet_CopanoBay, 0, b3

BullShark_MesquiteBay -> AllMenhaden_MesquiteBay, 0, b4
BullShark_AransasBay -> AllMenhaden_AransasBay, 0, b5
BullShark_CopanoBay -> AllMenhaden_CopanoBay, 0, b6

AlligatorGar_MesquiteBay -> AllMullet_MesquiteBay, 0, b7
AlligatorGar_AransasBay -> AllMullet_AransasBay, 0, b8
AlligatorGar_CopanoBay -> AllMullet_CopanoBay, 0, b9

AlligatorGar_MesquiteBay -> AllMenhaden_MesquiteBay, 0, b10
AlligatorGar_AransasBay -> AllMenhaden_AransasBay, 0, b11
AlligatorGar_CopanoBay -> AllMenhaden_CopanoBay, 0, b12

# Lagged relationships for predators (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888

BullShark_MesquiteBay -> BullShark_MesquiteBay, 1, l25
BullShark_AransasBay -> BullShark_AransasBay, 1, l26
BullShark_CopanoBay -> BullShark_CopanoBay, 1, l27

AlligatorGar_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l28
AlligatorGar_AransasBay -> AlligatorGar_AransasBay, 1, l29
AlligatorGar_CopanoBay -> AlligatorGar_CopanoBay, 1, l30

AllMullet_MesquiteBay -> AllMullet_MesquiteBay, 1, l13
AllMullet_AransasBay -> AllMullet_AransasBay, 1, l14
AllMullet_CopanoBay -> AllMullet_CopanoBay, 1, l15

AllMenhaden_MesquiteBay -> AllMenhaden_MesquiteBay, 1, l16
AllMenhaden_AransasBay -> AllMenhaden_AransasBay, 1, l17
AllMenhaden_CopanoBay -> AllMenhaden_CopanoBay, 1, l18
"

# Fit the model for prey prediction in Aransas Bay
fit_semAB_Pred_fulltopdown = dsem(sem = semAB_Pred_fulltopdown,
                      tsdata = AB_Pred_TS,
                      control = dsem_control(
                        quiet = TRUE))
summary(fit_semAB_Pred_fulltopdown)

# top down no density dependence model
semAB_Pred_topdown_noDD <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Predators predicting prey
BullShark_MesquiteBay -> AllMullet_MesquiteBay, 0, b1
BullShark_AransasBay -> AllMullet_AransasBay, 0, b2
BullShark_CopanoBay -> AllMullet_CopanoBay, 0, b3

BullShark_MesquiteBay -> AllMenhaden_MesquiteBay, 0, b4
BullShark_AransasBay -> AllMenhaden_AransasBay, 0, b5
BullShark_CopanoBay -> AllMenhaden_CopanoBay, 0, b6

AlligatorGar_MesquiteBay -> AllMullet_MesquiteBay, 0, b7
AlligatorGar_AransasBay -> AllMullet_AransasBay, 0, b8
AlligatorGar_CopanoBay -> AllMullet_CopanoBay, 0, b9

AlligatorGar_MesquiteBay -> AllMenhaden_MesquiteBay, 0, b10
AlligatorGar_AransasBay -> AllMenhaden_AransasBay, 0, b11
AlligatorGar_CopanoBay -> AllMenhaden_CopanoBay, 0, b12

# Lagged relationships for predators (Lag = 2)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888
"

# Fit the model for prey prediction in Aransas Bay
fit_semAB_Pred_topdown_noDD = dsem(sem = semAB_Pred_topdown_noDD ,
                           tsdata = AB_Pred_TS,
                           control = dsem_control(
                             quiet = TRUE))
summary(fit_semAB_Pred_topdown_noDD)

AIC(fit_semAB_Pred_abiotic)
AIC(fit_semAB_Pred_notrophics)
AIC(fit_semAB_Pred_topdown_noDD)
AIC(fit_semAB_Pred_bottomupnoDD)
AIC(fit_semAB_Pred_fulltopdown) #WINNER
AIC(fit_semAB_Pred_fullbottomup)

# Export model parameters (eg p values) as excel file
table_semAB_Pred<-summary(fit_semAB_Pred_fulltopdown)
write_xlsx(table_semAB_Pred, path = "table_semAB_Pred.xlsx")



###################################################################



# get data into a time series object from all numeric columns except YEAR
GB_Sciaenid_TS <- ts(
  GalvestonBay_Sciaenid_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(GalvestonBay_Sciaenid_Wide_Trans$YEAR)),  
  end = c(max(GalvestonBay_Sciaenid_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

# specify the dsem model for the sciaenid system for Galveston Bay for zero and 1 year lags

# full bottom up
semGB_Sciaenid_fullbottomup <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Atlanticcroaker_TrinityBay, 0, p5
PDSI -> Atlanticcroaker_GalvestonBay, 0, p6
PDSI -> Atlanticcroaker_WestBay, 0, p7
PDSI -> Atlanticcroaker_EastBay, 0, p8

PDSI -> BlueCrabSmall_TrinityBay, 0, p9
PDSI -> BlueCrabSmall_GalvestonBay, 0, p10
PDSI -> BlueCrabSmall_WestBay, 0, p11
PDSI -> BlueCrabSmall_EastBay, 0, p12

PDSI -> RedDrum_TrinityBay, 0, p13
PDSI -> RedDrum_GalvestonBay, 0, p14
PDSI -> RedDrum_WestBay, 0, p15
PDSI -> RedDrum_EastBay, 0, p16

PDSI -> SpottedSeatrout_TrinityBay, 0, p17
PDSI -> SpottedSeatrout_GalvestonBay, 0, p18
PDSI -> SpottedSeatrout_WestBay, 0, p19
PDSI -> SpottedSeatrout_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Atlanticcroaker_TrinityBay, 0, s1
Salinity_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, s2
Salinity_WestBay -> Atlanticcroaker_WestBay, 0, s3
Salinity_EastBay -> Atlanticcroaker_EastBay, 0, s4

Salinity_TrinityBay -> BlueCrabSmall_TrinityBay, 0, s5
Salinity_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, s6
Salinity_WestBay -> BlueCrabSmall_WestBay, 0, s7
Salinity_EastBay -> BlueCrabSmall_EastBay, 0, s8

Salinity_TrinityBay -> RedDrum_TrinityBay, 0, s9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 0, s10
Salinity_WestBay -> RedDrum_WestBay, 0, s11
Salinity_EastBay -> RedDrum_EastBay, 0, s12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 0, s13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, s14
Salinity_WestBay -> SpottedSeatrout_WestBay, 0, s15
Salinity_EastBay -> SpottedSeatrout_EastBay, 0, s16

# Prey predicting Sciaenid
Atlanticcroaker_TrinityBay -> RedDrum_TrinityBay, 0, b1
Atlanticcroaker_GalvestonBay -> RedDrum_GalvestonBay, 0, b2
Atlanticcroaker_WestBay -> RedDrum_WestBay, 0, b3
Atlanticcroaker_EastBay -> RedDrum_EastBay, 0, b4

BlueCrabSmall_TrinityBay -> RedDrum_TrinityBay, 0, b5
BlueCrabSmall_GalvestonBay -> RedDrum_GalvestonBay, 0, b6
BlueCrabSmall_WestBay -> RedDrum_WestBay, 0, b7
BlueCrabSmall_EastBay -> RedDrum_EastBay, 0, b8

Atlanticcroaker_TrinityBay -> SpottedSeatrout_TrinityBay, 0, b9
Atlanticcroaker_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, b10
Atlanticcroaker_WestBay -> SpottedSeatrout_WestBay, 0, b11
Atlanticcroaker_EastBay -> SpottedSeatrout_EastBay, 0, b12

BlueCrabSmall_TrinityBay -> SpottedSeatrout_TrinityBay, 0, b13
BlueCrabSmall_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, b14
BlueCrabSmall_WestBay -> SpottedSeatrout_WestBay, 0, b15
BlueCrabSmall_EastBay -> SpottedSeatrout_EastBay, 0, b16

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_TrinityBay, 1, l1
PDSI -> RedDrum_GalvestonBay, 1, l2
PDSI -> RedDrum_WestBay, 1, l3
PDSI -> RedDrum_EastBay, 1, l4

PDSI -> SpottedSeatrout_TrinityBay, 1, l5
PDSI -> SpottedSeatrout_GalvestonBay, 1, l6
PDSI -> SpottedSeatrout_WestBay, 1, l7
PDSI -> SpottedSeatrout_EastBay, 1, l8

RedDrum_TrinityBay -> RedDrum_TrinityBay, 1, l9
RedDrum_GalvestonBay -> RedDrum_GalvestonBay, 1, l10
RedDrum_WestBay -> RedDrum_WestBay, 1, l11
RedDrum_EastBay -> RedDrum_EastBay, 1, l12

SpottedSeatrout_TrinityBay -> SpottedSeatrout_TrinityBay, 1, l13
SpottedSeatrout_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, l14
SpottedSeatrout_WestBay -> SpottedSeatrout_WestBay, 1, l15
SpottedSeatrout_EastBay -> SpottedSeatrout_EastBay, 1, l16

Salinity_TrinityBay -> RedDrum_TrinityBay, 1, l17
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 1, l18
Salinity_WestBay -> RedDrum_WestBay, 1, l19
Salinity_EastBay -> RedDrum_EastBay, 1, l20

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 1, l21
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, l22
Salinity_WestBay -> SpottedSeatrout_WestBay, 1, l23
Salinity_EastBay -> SpottedSeatrout_EastBay, 1, l24
"

fit_semGB_Sciaenid_fullbottomup = dsem(sem = semGB_Sciaenid_fullbottomup,
                                       tsdata = GB_Sciaenid_TS,
                                       control = dsem_control(
                                         quiet = TRUE))
summary(fit_semGB_Sciaenid_fullbottomup)

# no DD bottom up model
semGB_Sciaenid_bottomup_noDD <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Atlanticcroaker_TrinityBay, 0, p5
PDSI -> Atlanticcroaker_GalvestonBay, 0, p6
PDSI -> Atlanticcroaker_WestBay, 0, p7
PDSI -> Atlanticcroaker_EastBay, 0, p8

PDSI -> BlueCrabSmall_TrinityBay, 0, p9
PDSI -> BlueCrabSmall_GalvestonBay, 0, p10
PDSI -> BlueCrabSmall_WestBay, 0, p11
PDSI -> BlueCrabSmall_EastBay, 0, p12

PDSI -> RedDrum_TrinityBay, 0, p13
PDSI -> RedDrum_GalvestonBay, 0, p14
PDSI -> RedDrum_WestBay, 0, p15
PDSI -> RedDrum_EastBay, 0, p16

PDSI -> SpottedSeatrout_TrinityBay, 0, p17
PDSI -> SpottedSeatrout_GalvestonBay, 0, p18
PDSI -> SpottedSeatrout_WestBay, 0, p19
PDSI -> SpottedSeatrout_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Atlanticcroaker_TrinityBay, 0, s1
Salinity_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, s2
Salinity_WestBay -> Atlanticcroaker_WestBay, 0, s3
Salinity_EastBay -> Atlanticcroaker_EastBay, 0, s4

Salinity_TrinityBay -> BlueCrabSmall_TrinityBay, 0, s5
Salinity_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, s6
Salinity_WestBay -> BlueCrabSmall_WestBay, 0, s7
Salinity_EastBay -> BlueCrabSmall_EastBay, 0, s8

Salinity_TrinityBay -> RedDrum_TrinityBay, 0, s9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 0, s10
Salinity_WestBay -> RedDrum_WestBay, 0, s11
Salinity_EastBay -> RedDrum_EastBay, 0, s12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 0, s13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, s14
Salinity_WestBay -> SpottedSeatrout_WestBay, 0, s15
Salinity_EastBay -> SpottedSeatrout_EastBay, 0, s16

# Prey predicting Sciaenid
Atlanticcroaker_TrinityBay -> RedDrum_TrinityBay, 0, b1
Atlanticcroaker_GalvestonBay -> RedDrum_GalvestonBay, 0, b2
Atlanticcroaker_WestBay -> RedDrum_WestBay, 0, b3
Atlanticcroaker_EastBay -> RedDrum_EastBay, 0, b4

BlueCrabSmall_TrinityBay -> RedDrum_TrinityBay, 0, b5
BlueCrabSmall_GalvestonBay -> RedDrum_GalvestonBay, 0, b6
BlueCrabSmall_WestBay -> RedDrum_WestBay, 0, b7
BlueCrabSmall_EastBay -> RedDrum_EastBay, 0, b8

Atlanticcroaker_TrinityBay -> SpottedSeatrout_TrinityBay, 0, b9
Atlanticcroaker_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, b10
Atlanticcroaker_WestBay -> SpottedSeatrout_WestBay, 0, b11
Atlanticcroaker_EastBay -> SpottedSeatrout_EastBay, 0, b12

BlueCrabSmall_TrinityBay -> SpottedSeatrout_TrinityBay, 0, b13
BlueCrabSmall_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, b14
BlueCrabSmall_WestBay -> SpottedSeatrout_WestBay, 0, b15
BlueCrabSmall_EastBay -> SpottedSeatrout_EastBay, 0, b16

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_TrinityBay, 1, l1
PDSI -> RedDrum_GalvestonBay, 1, l2
PDSI -> RedDrum_WestBay, 1, l3
PDSI -> RedDrum_EastBay, 1, l4

PDSI -> SpottedSeatrout_TrinityBay, 1, l5
PDSI -> SpottedSeatrout_GalvestonBay, 1, l6
PDSI -> SpottedSeatrout_WestBay, 1, l7
PDSI -> SpottedSeatrout_EastBay, 1, l8

Salinity_TrinityBay -> RedDrum_TrinityBay, 1, l9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 1, l10
Salinity_WestBay -> RedDrum_WestBay, 1, l11
Salinity_EastBay -> RedDrum_EastBay, 1, l12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 1, l13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, l14
Salinity_WestBay -> SpottedSeatrout_WestBay, 1, l15
Salinity_EastBay -> SpottedSeatrout_EastBay, 1, l16
"

fit_semGB_Sciaenid_bottomup_noDD = dsem(sem = semGB_Sciaenid_bottomup_noDD,
                                        tsdata = GB_Sciaenid_TS,
                                        control = dsem_control(
                                          quiet = TRUE))

summary(fit_semGB_Sciaenid_bottomup_noDD)

# abiotic only model
semGB_Sciaenid_abiotic <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Atlanticcroaker_TrinityBay, 0, p5
PDSI -> Atlanticcroaker_GalvestonBay, 0, p6
PDSI -> Atlanticcroaker_WestBay, 0, p7
PDSI -> Atlanticcroaker_EastBay, 0, p8

PDSI -> BlueCrabSmall_TrinityBay, 0, p9
PDSI -> BlueCrabSmall_GalvestonBay, 0, p10
PDSI -> BlueCrabSmall_WestBay, 0, p11
PDSI -> BlueCrabSmall_EastBay, 0, p12

PDSI -> RedDrum_TrinityBay, 0, p13
PDSI -> RedDrum_GalvestonBay, 0, p14
PDSI -> RedDrum_WestBay, 0, p15
PDSI -> RedDrum_EastBay, 0, p16

PDSI -> SpottedSeatrout_TrinityBay, 0, p17
PDSI -> SpottedSeatrout_GalvestonBay, 0, p18
PDSI -> SpottedSeatrout_WestBay, 0, p19
PDSI -> SpottedSeatrout_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Atlanticcroaker_TrinityBay, 0, s1
Salinity_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, s2
Salinity_WestBay -> Atlanticcroaker_WestBay, 0, s3
Salinity_EastBay -> Atlanticcroaker_EastBay, 0, s4

Salinity_TrinityBay -> BlueCrabSmall_TrinityBay, 0, s5
Salinity_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, s6
Salinity_WestBay -> BlueCrabSmall_WestBay, 0, s7
Salinity_EastBay -> BlueCrabSmall_EastBay, 0, s8

Salinity_TrinityBay -> RedDrum_TrinityBay, 0, s9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 0, s10
Salinity_WestBay -> RedDrum_WestBay, 0, s11
Salinity_EastBay -> RedDrum_EastBay, 0, s12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 0, s13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, s14
Salinity_WestBay -> SpottedSeatrout_WestBay, 0, s15
Salinity_EastBay -> SpottedSeatrout_EastBay, 0, s16

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_TrinityBay, 1, l1
PDSI -> RedDrum_GalvestonBay, 1, l2
PDSI -> RedDrum_WestBay, 1, l3
PDSI -> RedDrum_EastBay, 1, l4

PDSI -> SpottedSeatrout_TrinityBay, 1, l5
PDSI -> SpottedSeatrout_GalvestonBay, 1, l6
PDSI -> SpottedSeatrout_WestBay, 1, l7
PDSI -> SpottedSeatrout_EastBay, 1, l8

Salinity_TrinityBay -> RedDrum_TrinityBay, 1, l9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 1, l10
Salinity_WestBay -> RedDrum_WestBay, 1, l11
Salinity_EastBay -> RedDrum_EastBay, 1, l12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 1, l13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, l14
Salinity_WestBay -> SpottedSeatrout_WestBay, 1, l15
Salinity_EastBay -> SpottedSeatrout_EastBay, 1, l16
"

fit_semGB_Sciaenid_abiotic = dsem(sem = semGB_Sciaenid_abiotic,
                                  tsdata = GB_Sciaenid_TS,
                                  control = dsem_control(
                                    quiet = TRUE))

summary(fit_semGB_Sciaenid_abiotic)

# Full top down model
semGB_Sciaenid_fulltopdown <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Atlanticcroaker_TrinityBay, 0, p5
PDSI -> Atlanticcroaker_GalvestonBay, 0, p6
PDSI -> Atlanticcroaker_WestBay, 0, p7
PDSI -> Atlanticcroaker_EastBay, 0, p8

PDSI -> BlueCrabSmall_TrinityBay, 0, p9
PDSI -> BlueCrabSmall_GalvestonBay, 0, p10
PDSI -> BlueCrabSmall_WestBay, 0, p11
PDSI -> BlueCrabSmall_EastBay, 0, p12

PDSI -> RedDrum_TrinityBay, 0, p13
PDSI -> RedDrum_GalvestonBay, 0, p14
PDSI -> RedDrum_WestBay, 0, p15
PDSI -> RedDrum_EastBay, 0, p16

PDSI -> SpottedSeatrout_TrinityBay, 0, p17
PDSI -> SpottedSeatrout_GalvestonBay, 0, p18
PDSI -> SpottedSeatrout_WestBay, 0, p19
PDSI -> SpottedSeatrout_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Atlanticcroaker_TrinityBay, 0, s1
Salinity_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, s2
Salinity_WestBay -> Atlanticcroaker_WestBay, 0, s3
Salinity_EastBay -> Atlanticcroaker_EastBay, 0, s4

Salinity_TrinityBay -> BlueCrabSmall_TrinityBay, 0, s5
Salinity_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, s6
Salinity_WestBay -> BlueCrabSmall_WestBay, 0, s7
Salinity_EastBay -> BlueCrabSmall_EastBay, 0, s8

Salinity_TrinityBay -> RedDrum_TrinityBay, 0, s9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 0, s10
Salinity_WestBay -> RedDrum_WestBay, 0, s11
Salinity_EastBay -> RedDrum_EastBay, 0, s12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 0, s13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, s14
Salinity_WestBay -> SpottedSeatrout_WestBay, 0, s15
Salinity_EastBay -> SpottedSeatrout_EastBay, 0, s16

# Sciaenid predicting prey
RedDrum_TrinityBay -> Atlanticcroaker_TrinityBay, 0, b1
RedDrum_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, b2
RedDrum_WestBay -> Atlanticcroaker_WestBay, 0, b3
RedDrum_EastBay -> Atlanticcroaker_EastBay, 0, b4

RedDrum_TrinityBay -> BlueCrabSmall_TrinityBay, 0, b5
RedDrum_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, b6
RedDrum_WestBay -> BlueCrabSmall_WestBay, 0, b7
RedDrum_EastBay -> BlueCrabSmall_EastBay, 0, b8

SpottedSeatrout_TrinityBay -> Atlanticcroaker_TrinityBay, 0, b9
SpottedSeatrout_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, b10
SpottedSeatrout_WestBay -> Atlanticcroaker_WestBay, 0, b11
SpottedSeatrout_EastBay -> Atlanticcroaker_EastBay, 0, b12

SpottedSeatrout_TrinityBay -> BlueCrabSmall_TrinityBay, 0, b13
SpottedSeatrout_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, b14
SpottedSeatrout_WestBay -> BlueCrabSmall_WestBay, 0, b15
SpottedSeatrout_EastBay -> BlueCrabSmall_EastBay, 0, b16

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_TrinityBay, 1, l1
PDSI -> RedDrum_GalvestonBay, 1, l2
PDSI -> RedDrum_WestBay, 1, l3
PDSI -> RedDrum_EastBay, 1, l4

PDSI -> SpottedSeatrout_TrinityBay, 1, l5
PDSI -> SpottedSeatrout_GalvestonBay, 1, l6
PDSI -> SpottedSeatrout_WestBay, 1, l7
PDSI -> SpottedSeatrout_EastBay, 1, l8

Salinity_TrinityBay -> RedDrum_TrinityBay, 1, l9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 1, l10
Salinity_WestBay -> RedDrum_WestBay, 1, l11
Salinity_EastBay -> RedDrum_EastBay, 1, l12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 1, l13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, l14
Salinity_WestBay -> SpottedSeatrout_WestBay, 1, l15
Salinity_EastBay -> SpottedSeatrout_EastBay, 1, l16

RedDrum_TrinityBay -> RedDrum_TrinityBay, 1, l17
RedDrum_GalvestonBay -> RedDrum_GalvestonBay, 1, l18
RedDrum_WestBay -> RedDrum_WestBay, 1, l19
RedDrum_EastBay -> RedDrum_EastBay, 1, l20

SpottedSeatrout_TrinityBay -> SpottedSeatrout_TrinityBay, 1, l21
SpottedSeatrout_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, l22
SpottedSeatrout_WestBay -> SpottedSeatrout_WestBay, 1, l23
SpottedSeatrout_EastBay -> SpottedSeatrout_EastBay, 1, l24
"

# Fit the model for prey prediction in Galveston Bay
fit_semGB_Sciaenid_fulltopdown = dsem(sem = semGB_Sciaenid_fulltopdown,
                                      tsdata = GB_Sciaenid_TS,
                                      control = dsem_control(
                                        quiet = TRUE))

summary(fit_semGB_Sciaenid_fulltopdown)

# Top down no DD model
semGB_Sciaenid_topdown_noDD <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Atlanticcroaker_TrinityBay, 0, p5
PDSI -> Atlanticcroaker_GalvestonBay, 0, p6
PDSI -> Atlanticcroaker_WestBay, 0, p7
PDSI -> Atlanticcroaker_EastBay, 0, p8

PDSI -> BlueCrabSmall_TrinityBay, 0, p9
PDSI -> BlueCrabSmall_GalvestonBay, 0, p10
PDSI -> BlueCrabSmall_WestBay, 0, p11
PDSI -> BlueCrabSmall_EastBay, 0, p12

PDSI -> RedDrum_TrinityBay, 0, p13
PDSI -> RedDrum_GalvestonBay, 0, p14
PDSI -> RedDrum_WestBay, 0, p15
PDSI -> RedDrum_EastBay, 0, p16

PDSI -> SpottedSeatrout_TrinityBay, 0, p17
PDSI -> SpottedSeatrout_GalvestonBay, 0, p18
PDSI -> SpottedSeatrout_WestBay, 0, p19
PDSI -> SpottedSeatrout_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Atlanticcroaker_TrinityBay, 0, s1
Salinity_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, s2
Salinity_WestBay -> Atlanticcroaker_WestBay, 0, s3
Salinity_EastBay -> Atlanticcroaker_EastBay, 0, s4

Salinity_TrinityBay -> BlueCrabSmall_TrinityBay, 0, s5
Salinity_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, s6
Salinity_WestBay -> BlueCrabSmall_WestBay, 0, s7
Salinity_EastBay -> BlueCrabSmall_EastBay, 0, s8

Salinity_TrinityBay -> RedDrum_TrinityBay, 0, s9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 0, s10
Salinity_WestBay -> RedDrum_WestBay, 0, s11
Salinity_EastBay -> RedDrum_EastBay, 0, s12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 0, s13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, s14
Salinity_WestBay -> SpottedSeatrout_WestBay, 0, s15
Salinity_EastBay -> SpottedSeatrout_EastBay, 0, s16

# Sciaenid predicting prey
RedDrum_TrinityBay -> Atlanticcroaker_TrinityBay, 0, b1
RedDrum_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, b2
RedDrum_WestBay -> Atlanticcroaker_WestBay, 0, b3
RedDrum_EastBay -> Atlanticcroaker_EastBay, 0, b4

RedDrum_TrinityBay -> BlueCrabSmall_TrinityBay, 0, b5
RedDrum_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, b6
RedDrum_WestBay -> BlueCrabSmall_WestBay, 0, b7
RedDrum_EastBay -> BlueCrabSmall_EastBay, 0, b8

SpottedSeatrout_TrinityBay -> Atlanticcroaker_TrinityBay, 0, b9
SpottedSeatrout_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, b10
SpottedSeatrout_WestBay -> Atlanticcroaker_WestBay, 0, b11
SpottedSeatrout_EastBay -> Atlanticcroaker_EastBay, 0, b12

SpottedSeatrout_TrinityBay -> BlueCrabSmall_TrinityBay, 0, b13
SpottedSeatrout_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, b14
SpottedSeatrout_WestBay -> BlueCrabSmall_WestBay, 0, b15
SpottedSeatrout_EastBay -> BlueCrabSmall_EastBay, 0, b16

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_TrinityBay, 1, l1
PDSI -> RedDrum_GalvestonBay, 1, l2
PDSI -> RedDrum_WestBay, 1, l3
PDSI -> RedDrum_EastBay, 1, l4

PDSI -> SpottedSeatrout_TrinityBay, 1, l5
PDSI -> SpottedSeatrout_GalvestonBay, 1, l6
PDSI -> SpottedSeatrout_WestBay, 1, l7
PDSI -> SpottedSeatrout_EastBay, 1, l8

Salinity_TrinityBay -> RedDrum_TrinityBay, 1, l9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 1, l10
Salinity_WestBay -> RedDrum_WestBay, 1, l11
Salinity_EastBay -> RedDrum_EastBay, 1, l12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 1, l13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, l14
Salinity_WestBay -> SpottedSeatrout_WestBay, 1, l15
Salinity_EastBay -> SpottedSeatrout_EastBay, 1, l16
"

# Fit the model for prey prediction in Galveston Bay
fit_semGB_Sciaenid_topdown_noDD  = dsem(sem = semGB_Sciaenid_topdown_noDD,
                                        tsdata = GB_Sciaenid_TS,
                                        control = dsem_control(
                                          quiet = TRUE))

summary(fit_semGB_Sciaenid_topdown_noDD)

# No Prey Model
semGB_Sciaenid_notrophics <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Atlanticcroaker_TrinityBay, 0, p5
PDSI -> Atlanticcroaker_GalvestonBay, 0, p6
PDSI -> Atlanticcroaker_WestBay, 0, p7
PDSI -> Atlanticcroaker_EastBay, 0, p8

PDSI -> BlueCrabSmall_TrinityBay, 0, p9
PDSI -> BlueCrabSmall_GalvestonBay, 0, p10
PDSI -> BlueCrabSmall_WestBay, 0, p11
PDSI -> BlueCrabSmall_EastBay, 0, p12

PDSI -> RedDrum_TrinityBay, 0, p13
PDSI -> RedDrum_GalvestonBay, 0, p14
PDSI -> RedDrum_WestBay, 0, p15
PDSI -> RedDrum_EastBay, 0, p16

PDSI -> SpottedSeatrout_TrinityBay, 0, p17
PDSI -> SpottedSeatrout_GalvestonBay, 0, p18
PDSI -> SpottedSeatrout_WestBay, 0, p19
PDSI -> SpottedSeatrout_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Atlanticcroaker_TrinityBay, 0, s1
Salinity_GalvestonBay -> Atlanticcroaker_GalvestonBay, 0, s2
Salinity_WestBay -> Atlanticcroaker_WestBay, 0, s3
Salinity_EastBay -> Atlanticcroaker_EastBay, 0, s4

Salinity_TrinityBay -> BlueCrabSmall_TrinityBay, 0, s5
Salinity_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, s6
Salinity_WestBay -> BlueCrabSmall_WestBay, 0, s7
Salinity_EastBay -> BlueCrabSmall_EastBay, 0, s8

Salinity_TrinityBay -> RedDrum_TrinityBay, 0, s9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 0, s10
Salinity_WestBay -> RedDrum_WestBay, 0, s11
Salinity_EastBay -> RedDrum_EastBay, 0, s12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 0, s13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, s14
Salinity_WestBay -> SpottedSeatrout_WestBay, 0, s15
Salinity_EastBay -> SpottedSeatrout_EastBay, 0, s16

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_TrinityBay, 1, l1
PDSI -> RedDrum_GalvestonBay, 1, l2
PDSI -> RedDrum_WestBay, 1, l3
PDSI -> RedDrum_EastBay, 1, l4

PDSI -> SpottedSeatrout_TrinityBay, 1, l5
PDSI -> SpottedSeatrout_GalvestonBay, 1, l6
PDSI -> SpottedSeatrout_WestBay, 1, l7
PDSI -> SpottedSeatrout_EastBay, 1, l8

Salinity_TrinityBay -> RedDrum_TrinityBay, 1, l9
Salinity_GalvestonBay -> RedDrum_GalvestonBay, 1, l10
Salinity_WestBay -> RedDrum_WestBay, 1, l11
Salinity_EastBay -> RedDrum_EastBay, 1, l12

Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 1, l13
Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, l14
Salinity_WestBay -> SpottedSeatrout_WestBay, 1, l15
Salinity_EastBay -> SpottedSeatrout_EastBay, 1, l16

RedDrum_TrinityBay -> RedDrum_TrinityBay, 1, l25
RedDrum_GalvestonBay -> RedDrum_GalvestonBay, 1, l26
RedDrum_WestBay -> RedDrum_WestBay, 1, l27
RedDrum_EastBay -> RedDrum_EastBay, 1, l28

SpottedSeatrout_TrinityBay -> SpottedSeatrout_TrinityBay, 1, l29
SpottedSeatrout_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, l30
SpottedSeatrout_WestBay -> SpottedSeatrout_WestBay, 1, l31
SpottedSeatrout_EastBay -> SpottedSeatrout_EastBay, 1, l32
"

# Fit the No Prey model
fit_semGB_Sciaenid_notrophics = dsem(sem = semGB_Sciaenid_notrophics,
                                     tsdata = GB_Sciaenid_TS,
                                     control = dsem_control(
                                       quiet = TRUE))

summary(fit_semGB_Sciaenid_notrophics)

# Compare AICs of all six models
AIC(fit_semGB_Sciaenid_abiotic)
AIC(fit_semGB_Sciaenid_notrophics) 
AIC(fit_semGB_Sciaenid_bottomup_noDD)
AIC(fit_semGB_Sciaenid_topdown_noDD)
AIC(fit_semGB_Sciaenid_fullbottomup) 
AIC(fit_semGB_Sciaenid_fulltopdown)#WINNER

# Export model parameters (eg p values) as excel file
table_semGB_Sciaenid<-summary(fit_semGB_Sciaenid_fulltopdown)
write_xlsx(table_semGB_Sciaenid, path = "table_semGB_Sciaenid.xlsx")

# Second System and Major Bay: Keystone Predator Trophic System and Galveston Bay



#get data into a time series object from all numeric columns except YEAR
GB_Pred_TS <- ts(
  GalvestonBay_Pred_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(GalvestonBay_Pred_Wide_Trans$YEAR)),  
  end = c(max(GalvestonBay_Pred_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

# Full bottom-up model
semGB_Pred_fullbottomup <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Mullet_TrinityBay, 0, p5
PDSI -> Mullet_GalvestonBay, 0, p6
PDSI -> Mullet_WestBay, 0, p7
PDSI -> Mullet_EastBay, 0, p8

PDSI -> Menhaden_TrinityBay, 0, p9
PDSI -> Menhaden_GalvestonBay, 0, p10
PDSI -> Menhaden_WestBay, 0, p11
PDSI -> Menhaden_EastBay, 0, p12

PDSI -> BullShark_TrinityBay, 0, p13
PDSI -> BullShark_GalvestonBay, 0, p14
PDSI -> BullShark_WestBay, 0, p15
PDSI -> BullShark_EastBay, 0, p16

PDSI -> AlligatorGar_TrinityBay, 0, p17
PDSI -> AlligatorGar_GalvestonBay, 0, p18
PDSI -> AlligatorGar_WestBay, 0, p19
PDSI -> AlligatorGar_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Mullet_TrinityBay, 0, s1
Salinity_GalvestonBay -> Mullet_GalvestonBay, 0, s2
Salinity_WestBay -> Mullet_WestBay, 0, s3
Salinity_EastBay -> Mullet_EastBay, 0, s4

Salinity_TrinityBay -> Menhaden_TrinityBay, 0, s5
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 0, s6
Salinity_WestBay -> Menhaden_WestBay, 0, s7
Salinity_EastBay -> Menhaden_EastBay, 0, s8

Salinity_TrinityBay -> BullShark_TrinityBay, 0, s9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 0, s10
Salinity_WestBay -> BullShark_WestBay, 0, s11
Salinity_EastBay -> BullShark_EastBay, 0, s12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 0, s13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 0, s14
Salinity_WestBay -> AlligatorGar_WestBay, 0, s15
Salinity_EastBay -> AlligatorGar_EastBay, 0, s16

# Prey predicting predators
Mullet_TrinityBay -> BullShark_TrinityBay, 0, b1
Mullet_GalvestonBay -> BullShark_GalvestonBay, 0, b2
Mullet_WestBay -> BullShark_WestBay, 0, b3
Mullet_EastBay -> BullShark_EastBay, 0, b4

Menhaden_TrinityBay -> BullShark_TrinityBay, 0, b5
Menhaden_GalvestonBay -> BullShark_GalvestonBay, 0, b6
Menhaden_WestBay -> BullShark_WestBay, 0, b7
Menhaden_EastBay -> BullShark_EastBay, 0, b8

Mullet_TrinityBay -> AlligatorGar_TrinityBay, 0, b9
Mullet_GalvestonBay -> AlligatorGar_GalvestonBay, 0, b10
Mullet_WestBay -> AlligatorGar_WestBay, 0, b11
Mullet_EastBay -> AlligatorGar_EastBay, 0, b12

Menhaden_TrinityBay -> AlligatorGar_TrinityBay, 0, b13
Menhaden_GalvestonBay -> AlligatorGar_GalvestonBay, 0, b14
Menhaden_WestBay -> AlligatorGar_WestBay, 0, b15
Menhaden_EastBay -> AlligatorGar_EastBay, 0, b16

# Lagged relationships (Lag = 1)
PDSI -> BullShark_TrinityBay, 1, l1
PDSI -> BullShark_GalvestonBay, 1, l2
PDSI -> BullShark_WestBay, 1, l3
PDSI -> BullShark_EastBay, 1, l4

PDSI -> AlligatorGar_TrinityBay, 1, l5
PDSI -> AlligatorGar_GalvestonBay, 1, l6
PDSI -> AlligatorGar_WestBay, 1, l7
PDSI -> AlligatorGar_EastBay, 1, l8

PDSI -> Mullet_TrinityBay, 1, l17
PDSI -> Mullet_GalvestonBay, 1, l18
PDSI -> Mullet_WestBay, 1, l19
PDSI -> Mullet_EastBay, 1, l20

PDSI -> Menhaden_TrinityBay, 1, l21
PDSI -> Menhaden_GalvestonBay, 1, l22
PDSI -> Menhaden_WestBay, 1, l23
PDSI -> Menhaden_EastBay, 1, l24

Salinity_TrinityBay -> BullShark_TrinityBay, 1, l9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 1, l10
Salinity_WestBay -> BullShark_WestBay, 1, l11
Salinity_EastBay -> BullShark_EastBay, 1, l12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 1, l13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 1, l14
Salinity_WestBay -> AlligatorGar_WestBay, 1, l15
Salinity_EastBay -> AlligatorGar_EastBay, 1, l16

Salinity_TrinityBay -> Mullet_TrinityBay, 1, l25
Salinity_GalvestonBay -> Mullet_GalvestonBay, 1, l26
Salinity_WestBay -> Mullet_WestBay, 1, l27
Salinity_EastBay -> Mullet_EastBay, 1, l28

Salinity_TrinityBay -> Menhaden_TrinityBay, 1, l29
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 1, l30
Salinity_WestBay -> Menhaden_WestBay, 1, l31
Salinity_EastBay -> Menhaden_EastBay, 1, l32

BullShark_TrinityBay -> BullShark_TrinityBay, 1, l25
BullShark_GalvestonBay -> BullShark_GalvestonBay, 1, l26
BullShark_WestBay -> BullShark_WestBay, 1, l27
BullShark_EastBay -> BullShark_EastBay, 1, l28

AlligatorGar_TrinityBay -> AlligatorGar_TrinityBay, 1, l29
AlligatorGar_GalvestonBay -> AlligatorGar_GalvestonBay, 1, l30
AlligatorGar_WestBay -> AlligatorGar_WestBay, 1, l31
AlligatorGar_EastBay -> AlligatorGar_EastBay, 1, l32

Mullet_TrinityBay -> Mullet_TrinityBay, 1, z9
Mullet_GalvestonBay -> Mullet_GalvestonBay, 1, z10
Mullet_WestBay -> Mullet_WestBay, 1, z11
Mullet_EastBay -> Mullet_EastBay, 1, z12

Menhaden_TrinityBay -> Menhaden_TrinityBay, 1, z13
Menhaden_GalvestonBay -> Menhaden_GalvestonBay, 1, z14
Menhaden_WestBay -> Menhaden_WestBay, 1, z15
Menhaden_EastBay -> Menhaden_EastBay, 1, z16
"

# Fit the Full Bottom-Up Model for Galveston Bay
fit_semGB_Pred_fullbottomup = dsem(sem = semGB_Pred_fullbottomup,
                                   tsdata = GB_Pred_TS,
                                   control = dsem_control(
                                     quiet = TRUE))
summary(fit_semGB_Pred_fullbottomup)

# Abiotic Only Model
semGB_Pred_abiotic <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Mullet_TrinityBay, 0, p5
PDSI -> Mullet_GalvestonBay, 0, p6
PDSI -> Mullet_WestBay, 0, p7
PDSI -> Mullet_EastBay, 0, p8

PDSI -> Menhaden_TrinityBay, 0, p9
PDSI -> Menhaden_GalvestonBay, 0, p10
PDSI -> Menhaden_WestBay, 0, p11
PDSI -> Menhaden_EastBay, 0, p12

PDSI -> BullShark_TrinityBay, 0, p13
PDSI -> BullShark_GalvestonBay, 0, p14
PDSI -> BullShark_WestBay, 0, p15
PDSI -> BullShark_EastBay, 0, p16

PDSI -> AlligatorGar_TrinityBay, 0, p17
PDSI -> AlligatorGar_GalvestonBay, 0, p18
PDSI -> AlligatorGar_WestBay, 0, p19
PDSI -> AlligatorGar_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Mullet_TrinityBay, 0, s1
Salinity_GalvestonBay -> Mullet_GalvestonBay, 0, s2
Salinity_WestBay -> Mullet_WestBay, 0, s3
Salinity_EastBay -> Mullet_EastBay, 0, s4

Salinity_TrinityBay -> Menhaden_TrinityBay, 0, s5
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 0, s6
Salinity_WestBay -> Menhaden_WestBay, 0, s7
Salinity_EastBay -> Menhaden_EastBay, 0, s8

Salinity_TrinityBay -> BullShark_TrinityBay, 0, s9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 0, s10
Salinity_WestBay -> BullShark_WestBay, 0, s11
Salinity_EastBay -> BullShark_EastBay, 0, s12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 0, s13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 0, s14
Salinity_WestBay -> AlligatorGar_WestBay, 0, s15
Salinity_EastBay -> AlligatorGar_EastBay, 0, s16

# Lagged relationships for predators (Lag = 1)
PDSI -> BullShark_TrinityBay, 1, l1
PDSI -> BullShark_GalvestonBay, 1, l2
PDSI -> BullShark_WestBay, 1, l3
PDSI -> BullShark_EastBay, 1, l4

PDSI -> AlligatorGar_TrinityBay, 1, l5
PDSI -> AlligatorGar_GalvestonBay, 1, l6
PDSI -> AlligatorGar_WestBay, 1, l7
PDSI -> AlligatorGar_EastBay, 1, l8

PDSI -> Mullet_TrinityBay, 1, l9
PDSI -> Mullet_GalvestonBay, 1, l10
PDSI -> Mullet_WestBay, 1, l11
PDSI -> Mullet_EastBay, 1, l12

PDSI -> Menhaden_TrinityBay, 1, l13
PDSI -> Menhaden_GalvestonBay, 1, l14
PDSI -> Menhaden_WestBay, 1, l15
PDSI -> Menhaden_EastBay, 1, l16

Salinity_TrinityBay -> Mullet_TrinityBay, 1, l17
Salinity_GalvestonBay -> Mullet_GalvestonBay, 1, l18
Salinity_WestBay -> Mullet_WestBay, 1, l19
Salinity_EastBay -> Mullet_EastBay, 1, l20

Salinity_TrinityBay -> Menhaden_TrinityBay, 1, l21
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 1, l22
Salinity_WestBay -> Menhaden_WestBay, 1, l23
Salinity_EastBay -> Menhaden_EastBay, 1, l24

Salinity_TrinityBay -> BullShark_TrinityBay, 1, l25
Salinity_GalvestonBay -> BullShark_GalvestonBay, 1, l26
Salinity_WestBay -> BullShark_WestBay, 1, l27
Salinity_EastBay -> BullShark_EastBay, 1, l28

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 1, l29
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 1, l30
Salinity_WestBay -> AlligatorGar_WestBay, 1, l31
Salinity_EastBay -> AlligatorGar_EastBay, 1, l32
"

fit_semGB_Pred_abiotic = dsem(sem = semGB_Pred_abiotic,
                              tsdata = GB_Pred_TS,
                              control = dsem_control(
                                quiet = TRUE))
summary(fit_semGB_Pred_abiotic)

# No Density Dependence bottom up Model
semGB_Pred_bottomupnoDD <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Mullet_TrinityBay, 0, p5
PDSI -> Mullet_GalvestonBay, 0, p6
PDSI -> Mullet_WestBay, 0, p7
PDSI -> Mullet_EastBay, 0, p8

PDSI -> Menhaden_TrinityBay, 0, p9
PDSI -> Menhaden_GalvestonBay, 0, p10
PDSI -> Menhaden_WestBay, 0, p11
PDSI -> Menhaden_EastBay, 0, p12

PDSI -> BullShark_TrinityBay, 0, p13
PDSI -> BullShark_GalvestonBay, 0, p14
PDSI -> BullShark_WestBay, 0, p15
PDSI -> BullShark_EastBay, 0, p16

PDSI -> AlligatorGar_TrinityBay, 0, p17
PDSI -> AlligatorGar_GalvestonBay, 0, p18
PDSI -> AlligatorGar_WestBay, 0, p19
PDSI -> AlligatorGar_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Mullet_TrinityBay, 0, s1
Salinity_GalvestonBay -> Mullet_GalvestonBay, 0, s2
Salinity_WestBay -> Mullet_WestBay, 0, s3
Salinity_EastBay -> Mullet_EastBay, 0, s4

Salinity_TrinityBay -> Menhaden_TrinityBay, 0, s5
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 0, s6
Salinity_WestBay -> Menhaden_WestBay, 0, s7
Salinity_EastBay -> Menhaden_EastBay, 0, s8

Salinity_TrinityBay -> BullShark_TrinityBay, 0, s9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 0, s10
Salinity_WestBay -> BullShark_WestBay, 0, s11
Salinity_EastBay -> BullShark_EastBay, 0, s12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 0, s13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 0, s14
Salinity_WestBay -> AlligatorGar_WestBay, 0, s15
Salinity_EastBay -> AlligatorGar_EastBay, 0, s16

# Prey predicting predators
Mullet_TrinityBay -> BullShark_TrinityBay, 0, b1
Mullet_GalvestonBay -> BullShark_GalvestonBay, 0, b2
Mullet_WestBay -> BullShark_WestBay, 0, b3
Mullet_EastBay -> BullShark_EastBay, 0, b4

Menhaden_TrinityBay -> BullShark_TrinityBay, 0, b5
Menhaden_GalvestonBay -> BullShark_GalvestonBay, 0, b6
Menhaden_WestBay -> BullShark_WestBay, 0, b7
Menhaden_EastBay -> BullShark_EastBay, 0, b8

Mullet_TrinityBay -> AlligatorGar_TrinityBay, 0, b9
Mullet_GalvestonBay -> AlligatorGar_GalvestonBay, 0, b10
Mullet_WestBay -> AlligatorGar_WestBay, 0, b11
Mullet_EastBay -> AlligatorGar_EastBay, 0, b12

Menhaden_TrinityBay -> AlligatorGar_TrinityBay, 0, b13
Menhaden_GalvestonBay -> AlligatorGar_GalvestonBay, 0, b14
Menhaden_WestBay -> AlligatorGar_WestBay, 0, b15
Menhaden_EastBay -> AlligatorGar_EastBay, 0, b16

# Lagged relationships (Lag = 1)
PDSI -> BullShark_TrinityBay, 1, l1
PDSI -> BullShark_GalvestonBay, 1, l2
PDSI -> BullShark_WestBay, 1, l3
PDSI -> BullShark_EastBay, 1, l4

PDSI -> AlligatorGar_TrinityBay, 1, l5
PDSI -> AlligatorGar_GalvestonBay, 1, l6
PDSI -> AlligatorGar_WestBay, 1, l7
PDSI -> AlligatorGar_EastBay, 1, l8

PDSI -> Mullet_TrinityBay, 1, l17
PDSI -> Mullet_GalvestonBay, 1, l18
PDSI -> Mullet_WestBay, 1, l19
PDSI -> Mullet_EastBay, 1, l20

PDSI -> Menhaden_TrinityBay, 1, l21
PDSI -> Menhaden_GalvestonBay, 1, l22
PDSI -> Menhaden_WestBay, 1, l23
PDSI -> Menhaden_EastBay, 1, l24

Salinity_TrinityBay -> BullShark_TrinityBay, 1, l9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 1, l10
Salinity_WestBay -> BullShark_WestBay, 1, l11
Salinity_EastBay -> BullShark_EastBay, 1, l12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 1, l13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 1, l14
Salinity_WestBay -> AlligatorGar_WestBay, 1, l15
Salinity_EastBay -> AlligatorGar_EastBay, 1, l16

Salinity_TrinityBay -> Mullet_TrinityBay, 1, l25
Salinity_GalvestonBay -> Mullet_GalvestonBay, 1, l26
Salinity_WestBay -> Mullet_WestBay, 1, l27
Salinity_EastBay -> Mullet_EastBay, 1, l28

Salinity_TrinityBay -> Menhaden_TrinityBay, 1, l29
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 1, l30
Salinity_WestBay -> Menhaden_WestBay, 1, l31
Salinity_EastBay -> Menhaden_EastBay, 1, l32
"

# Fit the No Density Dependence bottom up model
fit_semGB_Pred_bottomupnoDD = dsem(sem = semGB_Pred_bottomupnoDD,
                                   tsdata = GB_Pred_TS,
                                   control = dsem_control(
                                     quiet = TRUE))
summary(fit_semGB_Pred_bottomupnoDD)


# No trophics Model
semGB_Pred_notrophics <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Mullet_TrinityBay, 0, p5
PDSI -> Mullet_GalvestonBay, 0, p6
PDSI -> Mullet_WestBay, 0, p7
PDSI -> Mullet_EastBay, 0, p8

PDSI -> Menhaden_TrinityBay, 0, p9
PDSI -> Menhaden_GalvestonBay, 0, p10
PDSI -> Menhaden_WestBay, 0, p11
PDSI -> Menhaden_EastBay, 0, p12

PDSI -> BullShark_TrinityBay, 0, p13
PDSI -> BullShark_GalvestonBay, 0, p14
PDSI -> BullShark_WestBay, 0, p15
PDSI -> BullShark_EastBay, 0, p16

PDSI -> AlligatorGar_TrinityBay, 0, p17
PDSI -> AlligatorGar_GalvestonBay, 0, p18
PDSI -> AlligatorGar_WestBay, 0, p19
PDSI -> AlligatorGar_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Mullet_TrinityBay, 0, s1
Salinity_GalvestonBay -> Mullet_GalvestonBay, 0, s2
Salinity_WestBay -> Mullet_WestBay, 0, s3
Salinity_EastBay -> Mullet_EastBay, 0, s4

Salinity_TrinityBay -> Menhaden_TrinityBay, 0, s5
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 0, s6
Salinity_WestBay -> Menhaden_WestBay, 0, s7
Salinity_EastBay -> Menhaden_EastBay, 0, s8

Salinity_TrinityBay -> BullShark_TrinityBay, 0, s9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 0, s10
Salinity_WestBay -> BullShark_WestBay, 0, s11
Salinity_EastBay -> BullShark_EastBay, 0, s12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 0, s13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 0, s14
Salinity_WestBay -> AlligatorGar_WestBay, 0, s15
Salinity_EastBay -> AlligatorGar_EastBay, 0, s16

# Lagged relationships (Lag = 1)
PDSI -> BullShark_TrinityBay, 1, l1
PDSI -> BullShark_GalvestonBay, 1, l2
PDSI -> BullShark_WestBay, 1, l3
PDSI -> BullShark_EastBay, 1, l4

PDSI -> AlligatorGar_TrinityBay, 1, l5
PDSI -> AlligatorGar_GalvestonBay, 1, l6
PDSI -> AlligatorGar_WestBay, 1, l7
PDSI -> AlligatorGar_EastBay, 1, l8

PDSI -> Mullet_TrinityBay, 1, l17
PDSI -> Mullet_GalvestonBay, 1, l18
PDSI -> Mullet_WestBay, 1, l19
PDSI -> Mullet_EastBay, 1, l20

PDSI -> Menhaden_TrinityBay, 1, l21
PDSI -> Menhaden_GalvestonBay, 1, l22
PDSI -> Menhaden_WestBay, 1, l23
PDSI -> Menhaden_EastBay, 1, l24

Salinity_TrinityBay -> BullShark_TrinityBay, 1, l9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 1, l10
Salinity_WestBay -> BullShark_WestBay, 1, l11
Salinity_EastBay -> BullShark_EastBay, 1, l12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 1, l13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 1, l14
Salinity_WestBay -> AlligatorGar_WestBay, 1, l15
Salinity_EastBay -> AlligatorGar_EastBay, 1, l16

Salinity_TrinityBay -> Mullet_TrinityBay, 1, l25
Salinity_GalvestonBay -> Mullet_GalvestonBay, 1, l26
Salinity_WestBay -> Mullet_WestBay, 1, l27
Salinity_EastBay -> Mullet_EastBay, 1, l28

Salinity_TrinityBay -> Menhaden_TrinityBay, 1, l29
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 1, l30
Salinity_WestBay -> Menhaden_WestBay, 1, l31
Salinity_EastBay -> Menhaden_EastBay, 1, l32

BullShark_TrinityBay -> BullShark_TrinityBay, 1, l25
BullShark_GalvestonBay -> BullShark_GalvestonBay, 1, l26
BullShark_WestBay -> BullShark_WestBay, 1, l27
BullShark_EastBay -> BullShark_EastBay, 1, l28

AlligatorGar_TrinityBay -> AlligatorGar_TrinityBay, 1, l29
AlligatorGar_GalvestonBay -> AlligatorGar_GalvestonBay, 1, l30
AlligatorGar_WestBay -> AlligatorGar_WestBay, 1, l31
AlligatorGar_EastBay -> AlligatorGar_EastBay, 1, l32

Mullet_TrinityBay -> Mullet_TrinityBay, 1, z9
Mullet_GalvestonBay -> Mullet_GalvestonBay, 1, z10
Mullet_WestBay -> Mullet_WestBay, 1, z11
Mullet_EastBay -> Mullet_EastBay, 1, z12

Menhaden_TrinityBay -> Menhaden_TrinityBay, 1, z13
Menhaden_GalvestonBay -> Menhaden_GalvestonBay, 1, z14
Menhaden_WestBay -> Menhaden_WestBay, 1, z15
Menhaden_EastBay -> Menhaden_EastBay, 1, z16
"

# Fit the No trophics SEM model
fit_semGB_Pred_notrophics = dsem(sem = semGB_Pred_notrophics,
                                 tsdata = GB_Pred_TS,
                                 control = dsem_control(
                                   quiet = TRUE))
summary(fit_semGB_Pred_notrophics)

# Full model for top down
semGB_Pred_fulltopdown <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Mullet_TrinityBay, 0, p5
PDSI -> Mullet_GalvestonBay, 0, p6
PDSI -> Mullet_WestBay, 0, p7
PDSI -> Mullet_EastBay, 0, p8

PDSI -> Menhaden_TrinityBay, 0, p9
PDSI -> Menhaden_GalvestonBay, 0, p10
PDSI -> Menhaden_WestBay, 0, p11
PDSI -> Menhaden_EastBay, 0, p12

PDSI -> BullShark_TrinityBay, 0, p13
PDSI -> BullShark_GalvestonBay, 0, p14
PDSI -> BullShark_WestBay, 0, p15
PDSI -> BullShark_EastBay, 0, p16

PDSI -> AlligatorGar_TrinityBay, 0, p17
PDSI -> AlligatorGar_GalvestonBay, 0, p18
PDSI -> AlligatorGar_WestBay, 0, p19
PDSI -> AlligatorGar_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Mullet_TrinityBay, 0, s1
Salinity_GalvestonBay -> Mullet_GalvestonBay, 0, s2
Salinity_WestBay -> Mullet_WestBay, 0, s3
Salinity_EastBay -> Mullet_EastBay, 0, s4

Salinity_TrinityBay -> Menhaden_TrinityBay, 0, s5
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 0, s6
Salinity_WestBay -> Menhaden_WestBay, 0, s7
Salinity_EastBay -> Menhaden_EastBay, 0, s8

Salinity_TrinityBay -> BullShark_TrinityBay, 0, s9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 0, s10
Salinity_WestBay -> BullShark_WestBay, 0, s11
Salinity_EastBay -> BullShark_EastBay, 0, s12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 0, s13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 0, s14
Salinity_WestBay -> AlligatorGar_WestBay, 0, s15
Salinity_EastBay -> AlligatorGar_EastBay, 0, s16

# Predators predicting prey
BullShark_TrinityBay -> Mullet_TrinityBay, 0, b1
BullShark_GalvestonBay -> Mullet_GalvestonBay, 0, b2
BullShark_WestBay -> Mullet_WestBay, 0, b3
BullShark_EastBay -> Mullet_EastBay, 0, b4

BullShark_TrinityBay -> Menhaden_TrinityBay, 0, b5
BullShark_GalvestonBay -> Menhaden_GalvestonBay, 0, b6
BullShark_WestBay -> Menhaden_WestBay, 0, b7
BullShark_EastBay -> Menhaden_EastBay, 0, b8

AlligatorGar_TrinityBay -> Mullet_TrinityBay, 0, b9
AlligatorGar_GalvestonBay -> Mullet_GalvestonBay, 0, b10
AlligatorGar_WestBay -> Mullet_WestBay, 0, b11
AlligatorGar_EastBay -> Mullet_EastBay, 0, b12

AlligatorGar_TrinityBay -> Menhaden_TrinityBay, 0, b13
AlligatorGar_GalvestonBay -> Menhaden_GalvestonBay, 0, b14
AlligatorGar_WestBay -> Menhaden_WestBay, 0, b15
AlligatorGar_EastBay -> Menhaden_EastBay, 0, b16

# Lagged relationships (Lag = 1)
PDSI -> BullShark_TrinityBay, 1, l1
PDSI -> BullShark_GalvestonBay, 1, l2
PDSI -> BullShark_WestBay, 1, l3
PDSI -> BullShark_EastBay, 1, l4

PDSI -> AlligatorGar_TrinityBay, 1, l5
PDSI -> AlligatorGar_GalvestonBay, 1, l6
PDSI -> AlligatorGar_WestBay, 1, l7
PDSI -> AlligatorGar_EastBay, 1, l8

PDSI -> Mullet_TrinityBay, 1, l17
PDSI -> Mullet_GalvestonBay, 1, l18
PDSI -> Mullet_WestBay, 1, l19
PDSI -> Mullet_EastBay, 1, l20

PDSI -> Menhaden_TrinityBay, 1, l21
PDSI -> Menhaden_GalvestonBay, 1, l22
PDSI -> Menhaden_WestBay, 1, l23
PDSI -> Menhaden_EastBay, 1, l24

Salinity_TrinityBay -> BullShark_TrinityBay, 1, l9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 1, l10
Salinity_WestBay -> BullShark_WestBay, 1, l11
Salinity_EastBay -> BullShark_EastBay, 1, l12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 1, l13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 1, l14
Salinity_WestBay -> AlligatorGar_WestBay, 1, l15
Salinity_EastBay -> AlligatorGar_EastBay, 1, l16

Salinity_TrinityBay -> Mullet_TrinityBay, 1, l25
Salinity_GalvestonBay -> Mullet_GalvestonBay, 1, l26
Salinity_WestBay -> Mullet_WestBay, 1, l27
Salinity_EastBay -> Mullet_EastBay, 1, l28

Salinity_TrinityBay -> Menhaden_TrinityBay, 1, l29
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 1, l30
Salinity_WestBay -> Menhaden_WestBay, 1, l31
Salinity_EastBay -> Menhaden_EastBay, 1, l32

BullShark_TrinityBay -> BullShark_TrinityBay, 1, l25
BullShark_GalvestonBay -> BullShark_GalvestonBay, 1, l26
BullShark_WestBay -> BullShark_WestBay, 1, l27
BullShark_EastBay -> BullShark_EastBay, 1, l28

AlligatorGar_TrinityBay -> AlligatorGar_TrinityBay, 1, l29
AlligatorGar_GalvestonBay -> AlligatorGar_GalvestonBay, 1, l30
AlligatorGar_WestBay -> AlligatorGar_WestBay, 1, l31
AlligatorGar_EastBay -> AlligatorGar_EastBay, 1, l32

Mullet_TrinityBay -> Mullet_TrinityBay, 1, z9
Mullet_GalvestonBay -> Mullet_GalvestonBay, 1, z10
Mullet_WestBay -> Mullet_WestBay, 1, z11
Mullet_EastBay -> Mullet_EastBay, 1, z12

Menhaden_TrinityBay -> Menhaden_TrinityBay, 1, z13
Menhaden_GalvestonBay -> Menhaden_GalvestonBay, 1, z14
Menhaden_WestBay -> Menhaden_WestBay, 1, z15
Menhaden_EastBay -> Menhaden_EastBay, 1, z16
"

# Fit the model for prey prediction in Galveston Bay
fit_semGB_Pred_fulltopdown = dsem(sem = semGB_Pred_fulltopdown,
                                  tsdata = GB_Pred_TS,
                                  control = dsem_control(
                                    quiet = TRUE))
summary(fit_semGB_Pred_fulltopdown)

# top down no density dependence model
semGB_Pred_topdown_noDD <- "
PDSI -> Salinity_TrinityBay, 0, p1
PDSI -> Salinity_GalvestonBay, 0, p2
PDSI -> Salinity_WestBay, 0, p3
PDSI -> Salinity_EastBay, 0, p4

PDSI -> Mullet_TrinityBay, 0, p5
PDSI -> Mullet_GalvestonBay, 0, p6
PDSI -> Mullet_WestBay, 0, p7
PDSI -> Mullet_EastBay, 0, p8

PDSI -> Menhaden_TrinityBay, 0, p9
PDSI -> Menhaden_GalvestonBay, 0, p10
PDSI -> Menhaden_WestBay, 0, p11
PDSI -> Menhaden_EastBay, 0, p12

PDSI -> BullShark_TrinityBay, 0, p13
PDSI -> BullShark_GalvestonBay, 0, p14
PDSI -> BullShark_WestBay, 0, p15
PDSI -> BullShark_EastBay, 0, p16

PDSI -> AlligatorGar_TrinityBay, 0, p17
PDSI -> AlligatorGar_GalvestonBay, 0, p18
PDSI -> AlligatorGar_WestBay, 0, p19
PDSI -> AlligatorGar_EastBay, 0, p20

# Salinity as a predictor
Salinity_TrinityBay -> Mullet_TrinityBay, 0, s1
Salinity_GalvestonBay -> Mullet_GalvestonBay, 0, s2
Salinity_WestBay -> Mullet_WestBay, 0, s3
Salinity_EastBay -> Mullet_EastBay, 0, s4

Salinity_TrinityBay -> Menhaden_TrinityBay, 0, s5
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 0, s6
Salinity_WestBay -> Menhaden_WestBay, 0, s7
Salinity_EastBay -> Menhaden_EastBay, 0, s8

Salinity_TrinityBay -> BullShark_TrinityBay, 0, s9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 0, s10
Salinity_WestBay -> BullShark_WestBay, 0, s11
Salinity_EastBay -> BullShark_EastBay, 0, s12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 0, s13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 0, s14
Salinity_WestBay -> AlligatorGar_WestBay, 0, s15
Salinity_EastBay -> AlligatorGar_EastBay, 0, s16

# Predators predicting prey
BullShark_TrinityBay -> Mullet_TrinityBay, 0, b1
BullShark_GalvestonBay -> Mullet_GalvestonBay, 0, b2
BullShark_WestBay -> Mullet_WestBay, 0, b3
BullShark_EastBay -> Mullet_EastBay, 0, b4

BullShark_TrinityBay -> Menhaden_TrinityBay, 0, b5
BullShark_GalvestonBay -> Menhaden_GalvestonBay, 0, b6
BullShark_WestBay -> Menhaden_WestBay, 0, b7
BullShark_EastBay -> Menhaden_EastBay, 0, b8

AlligatorGar_TrinityBay -> Mullet_TrinityBay, 0, b9
AlligatorGar_GalvestonBay -> Mullet_GalvestonBay, 0, b10
AlligatorGar_WestBay -> Mullet_WestBay, 0, b11
AlligatorGar_EastBay -> Mullet_EastBay, 0, b12

AlligatorGar_TrinityBay -> Menhaden_TrinityBay, 0, b13
AlligatorGar_GalvestonBay -> Menhaden_GalvestonBay, 0, b14
AlligatorGar_WestBay -> Menhaden_WestBay, 0, b15
AlligatorGar_EastBay -> Menhaden_EastBay, 0, b16

# Lagged relationships (Lag = 1)
PDSI -> BullShark_TrinityBay, 1, l1
PDSI -> BullShark_GalvestonBay, 1, l2
PDSI -> BullShark_WestBay, 1, l3
PDSI -> BullShark_EastBay, 1, l4

PDSI -> AlligatorGar_TrinityBay, 1, l5
PDSI -> AlligatorGar_GalvestonBay, 1, l6
PDSI -> AlligatorGar_WestBay, 1, l7
PDSI -> AlligatorGar_EastBay, 1, l8

PDSI -> Mullet_TrinityBay, 1, l17
PDSI -> Mullet_GalvestonBay, 1, l18
PDSI -> Mullet_WestBay, 1, l19
PDSI -> Mullet_EastBay, 1, l20

PDSI -> Menhaden_TrinityBay, 1, l21
PDSI -> Menhaden_GalvestonBay, 1, l22
PDSI -> Menhaden_WestBay, 1, l23
PDSI -> Menhaden_EastBay, 1, l24

Salinity_TrinityBay -> BullShark_TrinityBay, 1, l9
Salinity_GalvestonBay -> BullShark_GalvestonBay, 1, l10
Salinity_WestBay -> BullShark_WestBay, 1, l11
Salinity_EastBay -> BullShark_EastBay, 1, l12

Salinity_TrinityBay -> AlligatorGar_TrinityBay, 1, l13
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 1, l14
Salinity_WestBay -> AlligatorGar_WestBay, 1, l15
Salinity_EastBay -> AlligatorGar_EastBay, 1, l16

Salinity_TrinityBay -> Mullet_TrinityBay, 1, l25
Salinity_GalvestonBay -> Mullet_GalvestonBay, 1, l26
Salinity_WestBay -> Mullet_WestBay, 1, l27
Salinity_EastBay -> Mullet_EastBay, 1, l28

Salinity_TrinityBay -> Menhaden_TrinityBay, 1, l29
Salinity_GalvestonBay -> Menhaden_GalvestonBay, 1, l30
Salinity_WestBay -> Menhaden_WestBay, 1, l31
Salinity_EastBay -> Menhaden_EastBay, 1, l32

"

# Fit the model for prey prediction in Aransas Bay
fit_semGB_Pred_topdown_noDD = dsem(sem = semGB_Pred_topdown_noDD ,
                                   tsdata = GB_Pred_TS,
                                   control = dsem_control(
                                     quiet = TRUE))
summary(fit_semGB_Pred_topdown_noDD)

AIC(fit_semGB_Pred_abiotic)
AIC(fit_semGB_Pred_notrophics)
AIC(fit_semGB_Pred_topdown_noDD)
AIC(fit_semGB_Pred_bottomupnoDD)
AIC(fit_semGB_Pred_fulltopdown) 
AIC(fit_semGB_Pred_fullbottomup) #WINNER

# Export model parameters (eg p values) as excel file
table_semGB_Pred<-summary(fit_semGB_Pred_fullbottomup)
write_xlsx(table_semGB_Pred, path = "table_semGB_Pred.xlsx")


# The full top-down model is the winner for Aransas Bay for the Keystone Predator System
# The full bottom-up model is the winner for Galveston Bay for the Keystone Predator System
# The full top-down model is the winner for Galveston Bay for the sciaenid system
# And the no-trophic relationship model is the winner for Aransas Bay for the sciaenid system 


