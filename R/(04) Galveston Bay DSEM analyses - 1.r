###############################################################################################################################
##
##  (04) Galveston Bay DSEM analyses - 1 
##
###############################################################################################################################

######## Load required libraries
library( dsem )
library( ggplot2 )
library( ggpubr )
library( ggraph )
library( phylopath )
library( dplyr )
library( ggdag )
library( readxl )
library( profvis )
library( network )
library( ggnetwork )
library( igraph )
library( ggnetwork )
library( ggrepel )
library( gridExtra )

######## Load required functions
source( make.filename( "DSEM utilities.R", DIR$Functions ) )
source( make.filename( "plot utility functions.R", DIR$Functions ) )

######## Define subfolders to store outputs and figures
Subfolder_name <- " Galveston Bay - 1"
if ( !dir.exists( paste0( DIR$Output, Subfolder_name ) ) ) { dir.create( paste0( DIR$Output, Subfolder_name ) ) } 
OutputDir <- paste0( DIR$Output, Subfolder_name )
if ( !dir.exists( paste0( DIR$Figures, Subfolder_name ) ) ) { dir.create( paste0( DIR$Figures, Subfolder_name ) ) } 
FiguresDir <- paste0( DIR$Figures, Subfolder_name )

######## Define some settings
#### Define cushion space so that labels do not extend past plot bounds when plotting networks
cushion <- 0.1  

######## Load and process data
#### Get mean fish values per year
file_name <- make.filename( file = "GalvestonBay_GN_cpue.xlsx", path = DIR$Input )
GalvestonBay_GN_cpue <- read_excel( file_name )
nrow( GalvestonBay_GN_cpue ) ## 3555 records
GalvestonBay_GN_cpue2 <- na.omit( GalvestonBay_GN_cpue )
nrow( GalvestonBay_GN_cpue2 ) ## 3464 records
GalvestonBay_GN_cpue2$MONTH <- substr( GalvestonBay_GN_cpue2$DATE, 6, 7 )
GalvestonBay_GN_cpue2$YEAR <- substr( GalvestonBay_GN_cpue2$DATE, 1, 4 )
GalvestonBay_GN_YearMean <- GalvestonBay_GN_cpue2 %>%
	group_by( YEAR ) %>%
	summarise( across( everything(), mean, na.rm = TRUE ) )
GalvestonBay_GN_YearMean$YEAR <- as.numeric( GalvestonBay_GN_YearMean$YEAR )

#### Get mean blue crab values per year
file_name <- make.filename( file = "GalvestonBay_BS_bluecrabs.xlsx", path = DIR$Input )
GalvestonBay_BS_bluecrabs <- read_excel( file_name )
nrow( GalvestonBay_BS_bluecrabs ) ## 8665 records
GalvestonBay_BS_bluecrabs2 <- na.omit( GalvestonBay_BS_bluecrabs )
nrow( GalvestonBay_BS_bluecrabs2 ) ## 8665 records
GalvestonBay_BLUES_YearMean <- GalvestonBay_BS_bluecrabs2 %>%
	group_by( YEAR ) %>%
	summarise( across( everything(), mean, na.rm = TRUE ) )
GalvestonBay_BLUES_YearMean$YEAR <- as.numeric( GalvestonBay_BLUES_YearMean$YEAR )
columns_to_merge <- GalvestonBay_BLUES_YearMean %>%
	dplyr::select( YEAR, BlueCrabSmall, BlueCrabMedium, BlueCrabLarge )

#### Get mean brown shrimp values  per year
file_name <- make.filename( file = "Galveston_BS_shrimps.xlsx", path = DIR$Input )
Galveston_BS_shrimps <- read_excel( file_name )
nrow( Galveston_BS_shrimps ) ## 8665 records
Galveston_BS_shrimps$YEAR <- substr( Galveston_BS_shrimps$DATE, 1, 4 )
GalvestonBay_BS_shrimps2 <- na.omit( Galveston_BS_shrimps )
GalvestonBay_shrimps_YearMean <- GalvestonBay_BS_shrimps2 %>%
	group_by( YEAR ) %>%
	summarise( across( everything(), mean, na.rm = TRUE ) )
GalvestonBay_shrimps_YearMean$YEAR <- as.numeric( GalvestonBay_shrimps_YearMean$YEAR )
columns_to_merge_shrimp <- GalvestonBay_shrimps_YearMean %>%
	dplyr::select( YEAR, BrownShrimp )

######## Merge data for fishes with data for blue crabs based on YEAR 
GalvestonBay_GN_YearMean <- GalvestonBay_GN_YearMean %>%
	left_join( columns_to_merge, by = "YEAR" )

######## Merge data for fishes and blue crabs with data for brown shrimp based on YEAR 
GalvestonBay_GN_YearMean <- GalvestonBay_GN_YearMean %>%
	left_join( columns_to_merge_shrimp, by = "YEAR" )

######## Remove 1983 from the data because there are no shrimp records for that year
GalvestonBay_GN_YearMean <- GalvestonBay_GN_YearMean %>%
	filter( YEAR != 1983 )

######## Further process the dataset
GalvestonBay_GN_YearMean_log <- GalvestonBay_GN_YearMean %>%
	mutate( across( c( BlueCrabSmall, Salinity, Temp, BlackDrum, 
		RedDrum, SpottedSeatrout, HardheadCatfish, GafftopsailCatfish, 
		BrownShrimp ), log ) ) %>%  #### Log transformation
	mutate( across( c( BlueCrabSmall, Salinity, Temp, BlackDrum, 
		RedDrum, SpottedSeatrout, HardheadCatfish, GafftopsailCatfish, 
		BrownShrimp ), ~ ( . - mean(.) ) / sd(.) ) ) #### Normalize to mean zero and standard deviation of one

######## Set up the sem structure 
sem1 <- "

	## No time lag
    	BlueCrabSmall ->  RedDrum, 0, a0
	Temp ->  RedDrum, 0, b0
	Temp ->  BlueCrabSmall, 0, c0
	BlueCrabSmall ->  BlackDrum, 0, f0
	Temp ->  BlackDrum, 0, g0
	BlueCrabSmall -> SpottedSeatrout, 0, i0
	Temp ->  SpottedSeatrout, 0, j0
	BlueCrabSmall -> HardheadCatfish, 0, l0
	Temp ->  HardheadCatfish, 0, m0
	BlueCrabSmall -> GafftopsailCatfish, 0, o0
	Temp ->  GafftopsailCatfish, 0, p0
	BrownShrimp ->  RedDrum, 0, r0
	Temp ->  BrownShrimp, 0, s0
	BrownShrimp ->  BlackDrum, 0, u0
	BrownShrimp -> SpottedSeatrout, 0, v0
	BrownShrimp -> HardheadCatfish, 0, w0
	BrownShrimp -> GafftopsailCatfish, 0, x0
	Salinity ->  RedDrum, 0, bb0
    	Salinity ->  BlueCrabSmall, 0, cc0
    	Salinity ->  BlackDrum, 0, gg0
    	Salinity ->  SpottedSeatrout, 0, jj0
    	Salinity ->  HardheadCatfish, 0, mm0
    	Salinity ->  GafftopsailCatfish, 0, pp0
    	Salinity ->  BrownShrimp, 0, ss0

	#### Lag of 1 year
	BlueCrabSmall ->  RedDrum, 1, a1
	Temp ->  RedDrum, 1, b1
	Temp ->  BlueCrabSmall, 1, c1
	BlueCrabSmall -> BlueCrabSmall, 1, d1
	RedDrum -> RedDrum, 1, e1
	BlueCrabSmall ->  BlackDrum, 1, f1
	Temp ->  BlackDrum, 1, g1
	BlackDrum -> BlackDrum, 1, h1
	BlueCrabSmall -> SpottedSeatrout, 1, i1
	Temp ->  SpottedSeatrout, 1, j1
	SpottedSeatrout -> SpottedSeatrout, 1, k1
	BlueCrabSmall -> HardheadCatfish, 1, l1
	Temp ->  HardheadCatfish, 1, m1
	HardheadCatfish -> HardheadCatfish, 1, n1
	BlueCrabSmall -> GafftopsailCatfish, 1, o1
	Temp ->  GafftopsailCatfish, 1, p1
	GafftopsailCatfish-> GafftopsailCatfish, 1, q1
	BrownShrimp ->  RedDrum, 1, r1
	Temp ->  BrownShrimp, 1, s1
	BrownShrimp -> BrownShrimp, 1, t1
	BrownShrimp ->  BlackDrum, 1, u1
	BrownShrimp -> SpottedSeatrout, 1, v1
	BrownShrimp -> HardheadCatfish, 1, w1
	BrownShrimp -> GafftopsailCatfish, 1, x1
	Salinity ->  RedDrum, 1, bb1
	Salinity ->  BlueCrabSmall, 1, cc1
	Salinity ->  BlackDrum, 1, gg1
	Salinity ->  SpottedSeatrout, 1, jj1
	Salinity ->  HardheadCatfish, 1, mm1
	Salinity ->  GafftopsailCatfish, 1, pp1
	Salinity ->  BrownShrimp, 1, ss1

"

######## Get variables of interest in a time series 
data_GB <- ts( GalvestonBay_GN_YearMean_log[,( c( "BlueCrabSmall", "Temp", "Salinity", "RedDrum", 
	"SpottedSeatrout", "BlackDrum", "HardheadCatfish", "GafftopsailCatfish", "BrownShrimp" ) )], 
	start = c( min( GalvestonBay_GN_YearMean_log$YEAR ) ) ) 

######## Fit the DSEM
fit <- dsem( sem = sem1,
	tsdata = data_GB,
	estimate_delta0 = TRUE,
	control = dsem_control(
		quiet = TRUE,
		getsd = TRUE, newton_loops = 0 ) ) 
summary( fit )

######## Plot some results
#### Plot results for the relationships with no time lag
p <- plot( as_fitted_DAG( fit, lag = 0 ), text_size = 4 ) +
	expand_limits( x = c( -0.2, 2 ), y = c( -0.2, 0 ) )
p
pause( 3 ); graphics.off()

#### Plot results for the relationships with a time lag of 1 year
p <- plot( as_fitted_DAG( fit, lag = 1 ), text_size = 4 ) +
	expand_limits( x = c( -0.2, 2 ), y = c( -0.2, 0 ) )
p
pause( 3 ); graphics.off()

######## Get a model fit summary
fit_summary <- summary( fit )

######## Define subsets of the model fit summary
fit_summary_filtered_0 <- fit_summary %>%
	select( first, second, Estimate, lag ) %>%  
	filter( first != second, lag == "0" ) 
fit_summary_filtered_1 <- fit_summary %>%
	select( first, second, Estimate, lag ) %>%  
	filter( first != second, lag == "1" ) 

######## Get model fit summaries in the appropriate format for ggnetwork
g_0 <- graph_from_data_frame( d = fit_summary_filtered_0, directed = TRUE )
net_0 <- ggnetwork( g_0 )
g_1 <- graph_from_data_frame( d = fit_summary_filtered_1, directed = TRUE )
net_1 <- ggnetwork( g_1 )

######## Process the "net" model fit summary objects - We want to show significant negative relationships as red(ish), 
######## significant positive relationships as blue(ish) and show non significant relationships with grey arrows
net_0 <- net_0 %>%
	mutate( edge_color = case_when( Estimate < -0.2 ~ "lightcoral", Estimate > 0.2 ~ "steelblue3", TRUE ~ "grey88" ),
		xend = x + ( xend - x ) * 0.95, yend = y + ( yend - y ) * 0.95 )
net_1 <- net_1 %>%
	mutate( edge_color = case_when( Estimate < -0.2 ~ "lightcoral", Estimate > 0.2 ~ "steelblue3", TRUE ~ "grey88" ),
		xend = x + ( xend - x ) * 0.95, yend = y + ( yend - y ) * 0.95 )

######## Produce a figure showing the significant paths with no lag and a lag of 1 year side by side
p0 <- ggplot( net_0, aes( x = x, y = y, xend = xend, yend = yend ) ) +
	geom_edges( aes( size = abs( Estimate ), color = edge_color ), 
		arrow = arrow( length = unit( 9, "pt" ), type = "open" ), curvature = 0.1 ) +
	geom_nodelabel( aes( label = name ), fill = "seashell2" ) +
	scale_size_continuous( range = c( 0.5, 2 ) ) +  
	scale_color_identity() +  
	theme_void() + 
	theme( legend.position = "none",
		plot.margin = margin( 10, 10, 10, 10, "pt" ) ) + 
	expand_limits( x = c( min( net_0$x ) - cushion, max( net_0$x ) + cushion), 
		y = c( min( net_0$y ) - cushion, max( net_0$y ) + cushion ) )

p1 <- ggplot( net_1, aes( x = x, y = y, xend = xend, yend = yend ) ) +
	geom_edges( aes( size = abs( Estimate ), color = edge_color ), 
		arrow = arrow( length = unit( 9, "pt" ), type = "open" ), curvature = 0.1 ) +
	geom_nodelabel( aes( label = name ), fill = "seashell2" ) +
	scale_size_continuous( range = c( 0.5, 2 ) ) +  
	scale_color_identity() +  
	theme_void() + 
	theme( legend.position = "none",
		plot.margin = margin( 1, 10, 10, 10, "pt" ) ) + 
	expand_limits( x = c( min( net_1$x ) - cushion, max( net_1$x ) + cushion), 
		y = c( min( net_1$y ) - cushion, max( net_1$y ) + cushion ) )

margin <- theme( plot.margin = unit( c( 0, 0.5, 0, 0.5 ), "cm" ) )
grid.arrange( grobs = lapply( list( p0, p1 ), "+", margin ), ncol = 2 )
pause( 3 ); graphics.off()

#### Produce the same figure, this time using the "plot_dsem_network" function
p0 <- plot_dsem_ggnetwork( fit, lag = "0" )
p1 <- plot_dsem_ggnetwork( fit, lag = "1" )
margin <- theme( plot.margin = unit( c( 0, 0.5, 0, 0.5 ), "cm" ) )
grid.arrange( grobs = lapply( list( p0, p1 ), "+", margin ), ncol = 2 )
pause( 3 ); graphics.off()

######## Evaluate the DSEM using DHARMa residuals
#### Calculate quantile residuals using the predictive distribution from 
#### a jacknife (i.e., leave-one-out predictive distribution)
fittedPredictedResponse <- loo_residuals( fit, what = "loo" )





########################################################################################################################
OLD CODE THAT WILL NEED FIXING
y_ir_ts <- simulate ( fit,
	nsim = 100,
	resimulate_gmrf = TRUE )
y_ir <- array( NA, dim = c( nrow( y_ir_ts[[1]] ), ncol( y_ir_ts[[1]] ), 100 ), dimnames = NULL)
for ( i in 1 : 100 ) { y_ir[,,i] <- unlist( y_ir_ts[[i]] ) }

#### Generate DHARMa residuals for blue small crab and plot the histogram and QQ-plot of these DHARMa residuals
variable <- "BlueCrabSmall"

## Generate DHARMa residuals
dharmaRes <- DHARMa::createDHARMa( simulatedResponse = y_ir[,1,,drop = T], 
	observedResponse = GalvestonBay_GN_YearMean_log$BlueCrabSmall, 
	fittedPredictedResponse = fit$obj$report()$y_tj[,1] )

## Produce and save an histogram of DHARMa residuals 
val <- dharmaRes$scaledResiduals
val[val == 0] <- -0.01
val[val == 1] <- 1.01
jpeg( make.filename( paste0( "Histogram_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		hist( val, breaks = seq( -0.02, 1.02, len = 23 ), col = c( "red", rep( "lightgrey", 20 ), "red" ),
			main = "", xlab = "Residuals", cex.main = 2.5 )
dev.off()

## Produce and save a QQ-plot of DHARMa residuals 
jpeg( make.filename( paste0( "QQplot_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		gap::qqunif( dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6,
			main = "", cex.main = 2.5 )
dev.off()

#### Generate DHARMa residuals for temperature and plot the histogram and QQ-plot of these DHARMa residuals
variable <- "Temp"

## Generate DHARMa residuals
dharmaRes <- DHARMa::createDHARMa( simulatedResponse = y_ir[,2,,drop = T], 
	observedResponse = GalvestonBay_GN_YearMean_log$Temp, 
	fittedPredictedResponse = fit$obj$report()$y_tj[,2] )

## Produce and save an histogram of DHARMa residuals 
val <- dharmaRes$scaledResiduals
val[val == 0] <- -0.01
val[val == 1] <- 1.01
jpeg( make.filename( paste0( "Histogram_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		hist( val, breaks = seq( -0.02, 1.02, len = 23 ), col = c( "red", rep( "lightgrey", 20 ), "red" ),
			main = "", xlab = "Residuals", cex.main = 2.5 )
dev.off()

## Produce and save a QQ-plot of DHARMa residuals 
jpeg( make.filename( paste0( "QQplot_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		gap::qqunif( dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6,
			main = "", cex.main = 2.5 )
dev.off()

#### Generate DHARMa residuals for salinity and plot the histogram and QQ-plot of these DHARMa residuals
variable <- "Salinity"

## Generate DHARMa residuals
dharmaRes <- DHARMa::createDHARMa( simulatedResponse = y_ir[,3,,drop = T], 
	observedResponse = GalvestonBay_GN_YearMean_log$Salinity, 
	fittedPredictedResponse = fit$obj$report()$y_tj[,3] )

## Produce and save an histogram of DHARMa residuals 
val <- dharmaRes$scaledResiduals
val[val == 0] <- -0.01
val[val == 1] <- 1.01
jpeg( make.filename( paste0( "Histogram_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		hist( val, breaks = seq( -0.02, 1.02, len = 23 ), col = c( "red", rep( "lightgrey", 20 ), "red" ),
			main = "", xlab = "Residuals", cex.main = 2.5 )
dev.off()

## Produce and save a QQ-plot of DHARMa residuals 
jpeg( make.filename( paste0( "QQplot_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		gap::qqunif( dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6,
			main = "", cex.main = 2.5 )
dev.off()

#### Generate DHARMa residuals for red drum and plot the histogram and QQ-plot of these DHARMa residuals
variable <- "RedDrum"

## Generate DHARMa residuals
dharmaRes <- DHARMa::createDHARMa( simulatedResponse = y_ir[,4,,drop = T], 
	observedResponse = GalvestonBay_GN_YearMean_log$RedDrum, 
	fittedPredictedResponse = fit$obj$report()$y_tj[,4] )

## Produce and save an histogram of DHARMa residuals 
val <- dharmaRes$scaledResiduals
val[val == 0] <- -0.01
val[val == 1] <- 1.01
jpeg( make.filename( paste0( "Histogram_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		hist( val, breaks = seq( -0.02, 1.02, len = 23 ), col = c( "red", rep( "lightgrey", 20 ), "red" ),
			main = "", xlab = "Residuals", cex.main = 2.5 )
dev.off()

## Produce and save a QQ-plot of DHARMa residuals 
jpeg( make.filename( paste0( "QQplot_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		gap::qqunif( dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6,
			main = "", cex.main = 2.5 )
dev.off()

#### Generate DHARMa residuals for spotted seatrout and plot the histogram and QQ-plot of these DHARMa residuals
variable <- "SpottedSeatrout"

## Generate DHARMa residuals
dharmaRes <- DHARMa::createDHARMa( simulatedResponse = y_ir[,5,,drop = T], 
	observedResponse = GalvestonBay_GN_YearMean_log$SpottedSeatrout, 
	fittedPredictedResponse = fit$obj$report()$y_tj[,5] )

## Produce and save an histogram of DHARMa residuals 
val <- dharmaRes$scaledResiduals
val[val == 0] <- -0.01
val[val == 1] <- 1.01
jpeg( make.filename( paste0( "Histogram_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		hist( val, breaks = seq( -0.02, 1.02, len = 23 ), col = c( "red", rep( "lightgrey", 20 ), "red" ),
			main = "", xlab = "Residuals", cex.main = 2.5 )
dev.off()

## Produce and save a QQ-plot of DHARMa residuals 
jpeg( make.filename( paste0( "QQplot_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		gap::qqunif( dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6,
			main = "", cex.main = 2.5 )
dev.off()

#### Generate DHARMa residuals for black drum and plot the histogram and QQ-plot of these DHARMa residuals
variable <- "BlackDrum"

## Generate DHARMa residuals
dharmaRes <- DHARMa::createDHARMa( simulatedResponse = y_ir[,6,,drop = T], 
	observedResponse = GalvestonBay_GN_YearMean_log$BlackDrum, 
	fittedPredictedResponse = fit$obj$report()$y_tj[,6] )

## Produce and save an histogram of DHARMa residuals 
val <- dharmaRes$scaledResiduals
val[val == 0] <- -0.01
val[val == 1] <- 1.01
jpeg( make.filename( paste0( "Histogram_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		hist( val, breaks = seq( -0.02, 1.02, len = 23 ), col = c( "red", rep( "lightgrey", 20 ), "red" ),
			main = "", xlab = "Residuals", cex.main = 2.5 )
dev.off()

## Produce and save a QQ-plot of DHARMa residuals 
jpeg( make.filename( paste0( "QQplot_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		gap::qqunif( dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6,
			main = "", cex.main = 2.5 )
dev.off()

#### Generate DHARMa residuals for hardhead catfish and plot the histogram and QQ-plot of these DHARMa residuals
variable <- "HardheadCatfish"

## Generate DHARMa residuals
dharmaRes <- DHARMa::createDHARMa( simulatedResponse = y_ir[,7,,drop = T], 
	observedResponse = GalvestonBay_GN_YearMean_log$HardheadCatfish, 
	fittedPredictedResponse = fit$obj$report()$y_tj[,7] )

## Produce and save an histogram of DHARMa residuals 
val <- dharmaRes$scaledResiduals
val[val == 0] <- -0.01
val[val == 1] <- 1.01
jpeg( make.filename( paste0( "Histogram_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		hist( val, breaks = seq( -0.02, 1.02, len = 23 ), col = c( "red", rep( "lightgrey", 20 ), "red" ),
			main = "", xlab = "Residuals", cex.main = 2.5 )
dev.off()

## Produce and save a QQ-plot of DHARMa residuals 
jpeg( make.filename( paste0( "QQplot_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		gap::qqunif( dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6,
			main = "", cex.main = 2.5 )
dev.off()

#### Generate DHARMa residuals for gafftopsail catfish and plot the histogram and QQ-plot of these DHARMa residuals
variable <- "GafftopsailCatfish"

## Generate DHARMa residuals
dharmaRes <- DHARMa::createDHARMa( simulatedResponse = y_ir[,8,,drop = T], 
	observedResponse = GalvestonBay_GN_YearMean_log$GafftopsailCatfish, 
	fittedPredictedResponse = fit$obj$report()$y_tj[,8] )

## Produce and save an histogram of DHARMa residuals 
val <- dharmaRes$scaledResiduals
val[val == 0] <- -0.01
val[val == 1] <- 1.01
jpeg( make.filename( paste0( "Histogram_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		hist( val, breaks = seq( -0.02, 1.02, len = 23 ), col = c( "red", rep( "lightgrey", 20 ), "red" ),
			main = "", xlab = "Residuals", cex.main = 2.5 )
dev.off()

## Produce and save a QQ-plot of DHARMa residuals 
jpeg( make.filename( paste0( "QQplot_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		gap::qqunif( dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6,
			main = "", cex.main = 2.5 )
dev.off()

#### Generate DHARMa residuals for brown shrimp and plot the histogram and QQ-plot of these DHARMa residuals
variable <- "BrownShrimp"

## Generate DHARMa residuals
dharmaRes <- DHARMa::createDHARMa( simulatedResponse = y_ir[,9,,drop = T], 
	observedResponse = GalvestonBay_GN_YearMean_log$BrownShrimp, 
	fittedPredictedResponse = fit$obj$report()$y_tj[,9] )

## Produce and save an histogram of DHARMa residuals 
val <- dharmaRes$scaledResiduals
val[val == 0] <- -0.01
val[val == 1] <- 1.01
jpeg( make.filename( paste0( "Histogram_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		hist( val, breaks = seq( -0.02, 1.02, len = 23 ), col = c( "red", rep( "lightgrey", 20 ), "red" ),
			main = "", xlab = "Residuals", cex.main = 2.5 )
dev.off()

## Produce and save a QQ-plot of DHARMa residuals 
jpeg( make.filename( paste0( "QQplot_of_DHARMa_residuals_", variable, ".png" ), FiguresDir ), 
	width = 6, height = 7, units = "in", res = 600 )
		gap::qqunif( dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6,
			main = "", cex.main = 2.5 )
dev.off()

