###############################################################################################################################
##
##  (02) Demonstration of dynamic factor analysis with dsem
##
###############################################################################################################################

######## Load required libraries
library( dsem )
library( MARSS )
library( ggplot2 )

######## Define helper function
#### Define the "grab" function
grab = \( x, name ) x[which( names( x ) == name )] 

#### Define a custom function to plot states
plot_states <- function( out, vars = 1 : ncol( out$tmb_inputs$data$y_tj ) ) {

	xhat_tj <- as.list( out$sdrep, report = TRUE, what = "Estimate" )$z_tj[,vars,drop = FALSE]
	xse_tj <- as.list( out$sdrep, report = TRUE, what = "Std. Error" )$z_tj[,vars,drop = FALSE]

	longform <- expand.grid( Year = time( tsdata ), Var = colnames( tsdata )[vars] )
	longform$est <- as.vector( xhat_tj )
	longform$se <- as.vector( xse_tj )
	longform$upper <- longform$est + 1.96 * longform$se
	longform$lower <- longform$est - 1.96 * longform$se
  	longform$data <- as.vector( tsdata[,vars,drop = FALSE] )
  
	ggplot( data = longform ) +  
    		geom_line( aes( x = Year, y = est ) ) +
		geom_point( aes( x = Year, y = data ), color = "blue", na.rm = TRUE ) +
		geom_ribbon( aes( ymax = as.numeric( upper ), ymin = as.numeric( lower ), x = Year ), 
			color = "grey", alpha = 0.2 ) + 
		facet_wrap( facets = vars( Var ), scales = "free", ncol = 2 )

}

######## Dynamic factor analysis with dsem
######## "dsem" is an R package for fitting dynamic structural equation models (DSEMs) with a simple 
######## user-interface and generic specification of simultaneous and lagged effects in a potentially 
######## recursive structure. Here, we highlight how DSEM can be used to implement dynamic factor analysis (DFA).  
######## We specifically replicate analysis using the Multivariate Autoregressive State-Space (MARSS) package, 
######## using data that are provided as an example in the MARSS package
data( harborSealWA, package = "MARSS" )

#### Define the number of factors
#### n_factors >= 3 does not seem to converge using DSEM or MARSS without penalties
n_factors <- 2 

#### Using MARSS
#### We first illustrate a DFA model using two factors, fitted using MARSS
## Process data
dat <- t( scale( harborSealWA[,c( "SJI", "EBays", "SJF", "PSnd", "HC" )] ) )

## DFA with 3 states; the BFGS method is used because it fits much faster for this model
fit_MARSS <- MARSS( dat, 
	 model = list( m = n_factors ),
	form = "dfa",
	method = "BFGS",
	silent = TRUE )

## We can then plot the estimated factors (latent variables)
## We plot states using all data
plot( fit_MARSS, plot.type = "xtT" )

## And the estimated predictor for measurements (manifest variables)
## We plot expectation for data using all data
plot( fit_MARSS, plot.type = "fitted.ytT" )

#### Full-rank covariance using DSEM 

#### In DSEM syntax, we can first fit a saturated (full-covariance) model using the argument "covs"
## Add factors to data
tsdata <- ts( cbind( harborSealWA[,c( "SJI", "EBays", "SJF", "PSnd", "HC" )] ), start = 1978 )

## Scale and center (matches with what MARSS does when fitting a DFA)
tsdata <- scale( tsdata, center = TRUE, scale = TRUE )

## Define the SEM
sem <- "
	## Random-walk process for variables 
  	SJF -> SJF, 1, NA, 1
  	SJI -> SJI, 1, NA, 1
  	EBays -> EBays, 1, NA, 1
  	PSnd -> PSnd, 1, NA, 1
  	HC -> HC, 1, NA, 1
"

## Initial fit
mydsem0 <- dsem( tsdata = tsdata,
	covs = c( "SJF, SJI, EBays, PSnd, HC" ),
	sem = sem,
	family = rep( "normal", 5 ),
	control = dsem_control( quiet = TRUE,
		run_model = FALSE ) )   

## Fix all measurement errors at diagonal and equal
map <- mydsem0$tmb_inputs$map
map$lnsigma_j <- factor( rep( 1, ncol( tsdata ) ) )

## Re-fit the model
mydsem_full <- dsem( tsdata = tsdata,
	covs = c( "SJF, SJI, EBays, PSnd, HC" ),
	sem = sem,
	family = rep( "normal", 5 ),
 	control = dsem_control( quiet = TRUE,
	map = map ) )

## Plot states
plot_states( mydsem_full )

## These estimated states follow the data more closely and have smaller estimated confidence intervals. 
## Presumably this occurs because we are using a full-rank covariance so far

#### Reduced-rank factor model with measurement errors

#### Next, we can specify two factors while eliminating additional process error and estimating measurement errors.  
#### This requires us to switch to gmrf_parameterization = "projection", so that we can fit a rank-deficient 
#### Gaussian Markov random field
## Add factors to data
tsdata <- harborSealWA[,c( "SJI", "EBays", "SJF", "PSnd", "HC" )]
newcols <- array( NA,
	dim = c( nrow( tsdata ), n_factors ),
	dimnames = list( NULL, paste0( "F", seq_len( n_factors ) ) ) )
tsdata <- ts( cbind( tsdata, newcols ), start = 1978 )

## Scale and center (matches with what MARSS does when fitting a DFA)
tsdata <- scale( tsdata, center = TRUE, scale = TRUE )

## Define the structural equation model
sem <- make_dfa( variables = c( "SJI", "EBays", "SJF", "PSnd", "HC" ),
	n_factors = n_factors )

## Initial fit
mydsem0 <- dsem( tsdata = tsdata,
	sem = sem,
	family = c( rep( "normal", 5 ), rep( "fixed", n_factors ) ),
	estimate_delta0 = TRUE,
	control = dsem_control( quiet = TRUE,
		run_model = FALSE,
		gmrf_parameterization = "projection" ) )

## Fix all measurement errors at diagonal and equal
map <- mydsem0$tmb_inputs$map
map$lnsigma_j <- factor( rep( 1, ncol( tsdata ) ) )

## Fix factors to have initial value, and variables to not
map$delta0_j <- factor( c( rep( NA, ncol( harborSealWA ) - 1 ), 1 : n_factors ) )

## Fix variables to have no stationary mean except what is predicted by initial value
map$mu_j <- factor( rep( NA, ncol( tsdata ) ) )

## Profile "delta0_j" to match MARSS (which treats initial condition as unpenalized random effect)
mydfa = dsem( tsdata = tsdata,
	sem = sem,
	family = c( rep( "normal", 5 ), rep( "fixed", n_factors ) ),
	estimate_delta0 = TRUE,
	control = dsem_control( quiet = TRUE,
		map = map,
		use_REML = TRUE,
		# profile = "delta0_j",
		gmrf_parameterization = "projection" ) )

## We again plot the estimated latent variables
## Plot estimated factors
plot_states( mydfa, vars = 5 + seq_len( n_factors ) )

## And the estimated predictor for manifest variables
## Plot estimated variables
plot_states( mydfa, vars = 1 : 5 )

## This results in similar (but not identical) factor values using MARSS and DSEM.  
## In particular, DSEM has higher variance in early years. This likely arises because the default MARSS implementation 
## of DFA includes a penalty of the initial state  with mean zero and variance. This term presumably provides additional 
## information about the initial year such that MARSS DFA results are not invariant to reversing the order of the data

#### To further explore, we can modify the MARSS DFA to eliminate the prior on initial conditions, 
#### based on help from Dr. Eli Holmes. This involves specifying the following
## Extract internal settings
modmats <- summary( fit_MARSS$model, silent = TRUE )

## Redefine defaults
modmats$V0 <- matrix( 0, n_factors, n_factors )
modmats$x0 <- "unequal" 

## Refit
fit_MARSS2 <- MARSS( dat, 
	model = modmats,
	silent = TRUE,
	control = list( abstol = 0.001,
		conv.test.slope.tol = 0.01, 
		maxit = 1000 ) )

## These have estimated time-series that are more similar to those from DSEM

## Plots states using all data
plot( fit_MARSS2, plot.type = "xtT" )

#### We can now compare the three options in terms of the fitted log-likelihood

#### Compare likelihood for MARSS and DSEM
Table <- c( "MARSS" = logLik( fit_MARSS ), 
	"DSEM" = logLik( mydfa ), 
	"MARSS_no_pen" = logLik( fit_MARSS2 ) )
knitr::kable( Table, digits = 3 )       

## This confirms that the MARSS model without a penalty on initial conditions results in the same likelihood as DSEM 

## Finally, we can also compare the three options in terms of estimated loadings
Table <- cbind( "MARSS" = as.vector( fit_MARSS$par$Z ),
	"DSEM" = grab( mydfa$opt$par, "beta_z" ),
	"MARSS_no_pen" = as.vector( fit_MARSS2$par$Z ) )
rownames( Table ) <- names( fit_MARSS$coef )[1 : nrow( Table )]
knitr::kable( Table, digits = 3 )       

## The estimating loadings are similar using DSEM and the MARSS model without initial penalty, 
## except with label switching (where some factors and loadings can be multiplied by -1 with no change in the model)
