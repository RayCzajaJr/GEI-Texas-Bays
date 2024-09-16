###############################################################################################################################
##
##  DSEM utility functions
##
###############################################################################################################################

######## Define the "plot_dsem_network" function
plot_dsem_network <- function ( fit, cushion = 0.1 ) {
  
	fit_summary <- summary( fit )
  
	fit_summary_filtered <- fit_summary %>%
		select( first, second, Estimate ) %>%  
		filter( first != second )
  
	g <- graph_from_data_frame( d = fit_summary_filtered, directed = TRUE )
  
	net <- ggnetwork( g )
  
	net <- net %>%
		mutate( edge_color = case_when( Estimate < -0.2 ~ "lightcoral", Estimate > 0.2 ~ "steelblue3",
			TRUE ~ "grey93" ),
      	xend = x + ( xend - x ) * 0.95,
      	yend = y + ( yend - y ) * 0.95 )
  
	ggplot( net, aes( x = x, y = y, xend = xend, yend = yend, layout = "kk" ) ) +
		geom_nodelabel( aes( label = name ), fill = "seashell2", size = 4 ) +
	geom_edges( aes( size = abs( Estimate ), color = edge_color, alpha = 0.9 ), 
		arrow = arrow( length = unit( 9, "pt" ), type = "open" ), curvature = 0.05 ) +
	scale_size_continuous( range = c( .0001, 2 ) ) + 
    	scale_color_identity() +  
    	theme_void() + 
    	theme( legend.position = "none", plot.margin = margin( 10, 10, 10, 10, "pt" ) ) +  
	expand_limits( x = c( min( net$x ) - cushion, max( net$x ) + cushion ), 
		y = c( min( net$y ) - cushion, max( net$y ) + cushion ) )

}

######## Define the "plot_dsem_ggnetwork" function
plot_dsem_ggnetwork <- function ( fit, lag, cushion = 0.1 ) {
  
	fit_summary <- summary( fit )
  
	fit_summary_filtered <- fit_summary %>%
		select( first, second, Estimate ) %>%  
		filter( first != second, lag == lag )
  
	g <- graph_from_data_frame( d = fit_summary_filtered, directed = TRUE )
  
	net <- ggnetwork( g )
  
	net <- net %>%
		mutate( edge_color = case_when( Estimate < -0.2 ~ "lightcoral", Estimate > 0.2 ~ "steelblue3",
			TRUE ~ "grey93" ),
      	xend = x + ( xend - x ) * 0.95,
      	yend = y + ( yend - y ) * 0.95 )
  
	ggplot( net, aes( x = x, y = y, xend = xend, yend = yend, layout = "kk" ) ) +
		geom_nodelabel( aes( label = name ), fill = "seashell2", size = 4 ) +
	geom_edges( aes( size = abs( Estimate ), color = edge_color, alpha = 0.9 ), 
		arrow = arrow( length = unit( 9, "pt" ), type = "open" ), curvature = 0.05 ) +
	scale_size_continuous( range = c( .0001, 2 ) ) + 
    	scale_color_identity() +  
    	theme_void() + 
    	theme( legend.position = "none", plot.margin = margin( 10, 10, 10, 10, "pt" ) ) +  
	expand_limits( x = c( min( net$x ) - cushion, max( net$x ) + cushion ), 
		y = c( min( net$y ) - cushion, max( net$y ) + cushion ) )

}