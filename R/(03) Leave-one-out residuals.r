###############################################################################################################################
##
##  (03) Demonstration of the calculation and use of leave-one-out residuals for dsem
##
###############################################################################################################################

######## Load required libraries
library( dsem )

######## Explore the calculation and use of leave-one-out residuals for dsem
x <- rnorm( 100 )
y <- 1 + 0.5 * x + rnorm( 100 )
x[sample( seq_along( x ), size = 10, replace = FALSE )] <- NA
y[sample( seq_along( x ), size = 10, replace = FALSE )] <- NA

tsdata <- ts( data.frame( x = x, y = y ) )

fit <- dsem( tsdata = tsdata,
            sem = "x -> y, 0, slope" )

resid <- loo_residuals( fit )

hist( resid )

#### Creates and plot a DHARMa object
samples <- loo_residuals( fit, what = "samples" ) 

which_use <- which(!is.na( tsdata ) )

fittedPredictedResponse <- loo_residuals( fit, what = "loo" )

res <- DHARMa::createDHARMa( simulatedResponse = apply( samples, MARGIN = 3, FUN = as.vector )[which_use,],
	observedResponse = as.vector( tsdata )[which_use],
	fittedPredictedResponse = fittedPredictedResponse$est )

plot( res )