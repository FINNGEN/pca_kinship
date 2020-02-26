aberrant <-
function ( x , lambda , niter = 10000 , burnin = 100 , prior_df = NULL, prior_scale= NULL, hyper_prior_mean = NULL , hyper_prior_var = NULL,  hyper_prior_df = NULL , hyper_prior_scale = NULL , alpha = NULL, beta = NULL , standardize = TRUE , verbose = TRUE , uncorr = FALSE )     
     {

     if ( missing ( x ) ) { stop ( "'x' missing" ) } 

     # n is the number of individuals, m is the number of summary statistics, and g is the number of groups (2 groups: "normal" individuals and outliers)
     x_ini <- as.matrix ( x ) 
     nb_ind <- nrow ( x_ini )
     m <- ncol ( x_ini )
     g <- 2

     ## Data and parameter format checks 
     # currently, the function can handle only 2 summary statistics
     if ( m != 2 ) 
	  {
	  stop ( "x must be a matrix of two summary statistics" )
	  }

     if ( missing ( lambda ) || ! is.numeric ( lambda ) || length ( lambda ) != 1 )
	  {
	  stop ("Parameter 'lambda' must be a numeric of length 1")
	  }

     if ( ! is.null ( hyper_prior_var ) & ( ! is.numeric ( hyper_prior_var ) || length ( hyper_prior_var  ) != m )  ) 
	  {
	  stop ( paste ( "Prior parameter 'hyper_prior_var' must be either NULL either a vector of length " , m ) ) 
	  }

     if ( ! is.null ( hyper_prior_df ) & ( ! is.numeric ( hyper_prior_df ) || length ( hyper_prior_df ) != 2 ) )
	  {
	  stop ( "Parameter 'hyper_prior_df' must be either NULL either a numeric vector of length 2" )
	  }

     if ( ! is.null ( hyper_prior_scale ) & ( ! is.numeric ( hyper_prior_scale ) || length ( hyper_prior_scale ) != m ) )
	  {
	  stop ( "Parameter 'hyper_prior_scale' must be either NULL either a numeric vector of length 2" )
	  }

     if ( ! is.null ( prior_df ) & ( ! is.numeric ( prior_df ) || ( uncorr == TRUE & length ( prior_df ) != 2 ) || ( uncorr == FALSE & length ( prior_df ) != 1 ) ) )
	  {
	  stop ( "Parameter 'prior_df' must be a numeric of length 1 if uncorr=FALSE and a numeric of length 2 if uncorr=TRUE" )
	  }

     if ( ! is.null ( prior_scale ) & ( ! is.matrix ( prior_scale ) || nrow ( prior_scale ) != m ) )
	  {
	  stop ( "Parameter 'prior_scale' must be either NULL either a 2*2 matrix" )
	  }

     if ( ! ( is.null ( hyper_prior_mean ) & is.null ( hyper_prior_var ) & is.null ( prior_df ) & is.null ( prior_scale ) & is.null ( hyper_prior_df ) & is.null ( hyper_prior_scale ) ) & ! ( ! is.null ( hyper_prior_mean ) &  ! is.null ( hyper_prior_var ) & ! is.null ( prior_df ) &  ! is.null ( prior_scale ) & !is.null ( hyper_prior_df ) & !is.null ( hyper_prior_scale ) ) )
	  {
	  stop ( paste ( "Prior parameters  'prior_df', 'prior_scale', 'hyper_prior_mean', 'hyper_prior_var', 'hyper_prior_df' and 'hyper_prior_scale' must be all NULL or all numeric" ) )
	  }
     
     if ( ! is.null (prior_scale) & uncorr == TRUE & ( ! prior_scale[1,2]==0 || ! prior_scale[2,1]==0 ) )
	  {
	  stop ( "If uncorr = TRUE, then matrix 'prior_scale' must be either NULL or diagonal" )
	  }

     if ( ! is.null ( alpha ) & ( ! is.numeric ( alpha ) || length ( alpha ) != 1 ) )
	  {
	  stop ( "Parameter 'alpha' must be a numeric of length 1" )
	  }

     if ( ! is.null ( beta ) & ( ! is.numeric ( beta ) || length ( beta ) != 1 ) )
	  {
	  stop ( "Parameter 'beta' must be a numeric of length 1" )
	  }

     if ( ! ( is.null ( alpha ) & is.null ( beta ) ) & ! ( ! is.null ( alpha ) &  ! is.null ( beta ) ) )
	  {
	  stop ( paste ( "Prior parameters 'alpha' and 'beta' must be both NULL or both numeric" ) )
	  }


     # Standardize the data if required
     if ( standardize ) { x <- scale ( x_ini ) }
     if ( ! standardize ) { x <- x_ini }

     # param is a list containing values for mean and variance of the summary statistics in the "normal" and "outlier" group. Hyper_param is a list containing values for mean and variance of the prior on the means of the distribution describing the variability of the summary statistics. 
     param <- list ( mean = array ( NA , dim = c ( m , g ) ) , var = array( 0 , dim = c ( m , m , g ) ) )
     hyper_param <- list ( mean = array ( NA , dim = m ) , var = array( 0 , dim = c ( m , m ) ) )

     # Initialization of various parameters
     # sum_mean: sum the sampled mean of the summary statistics for the "normal" individuals
     # sum_var: sum the sampled variance of the summary statistics for the "normal" individuals
     # group: 0 if individual is "normal", 1 otherwise
     # post: posterior probability of being an outlier for each individual
     # counts: number of iterations where each individual is sampled as an outliers	
     # theta: parameter of the Bernoulli distribution for the outlier indicator z

     sum_mean <- rep ( 0 , m )
     sum_var <- matrix ( 0 , m , m )
     group <- rep ( NA , nb_ind )
     post <- rep ( NA , nb_ind )
     counts <- rep ( 0 , nb_ind )
     theta <- rep ( 0 , nb_ind )

     ## Initialization of algorithm
     # if no prior is given, individuals are considered outliers if at least one of their summary statistics is higher (respectively lower) than the corresponding 90% quantile (respectively 10% quantile). Mean of outlier individuals is set to the mean of the "normal" individuals and covariance matrix of the outliers is set to lamda^2 times the covariance matrix of the "normal" individual
     # if some priors are given, then the parameters are set according to these priors (mean of both normal and outliers are set to hyper_prior_mean, and the variance of the normal individuals is set to the mode of the prior distribution used for this parameter)
     # out: index of "outliers"
     # n: number of "normal" and "outlier" individuals
     # q: prior probability of being an outlier
     if ( is.null ( alpha ) ) { q <- 1/2 }
     if ( ! is.null (alpha) ) { q <- alpha / (alpha + beta) }

     if ( is.null(hyper_prior_mean) )
	{
	out <- x [,1] <= quantile ( x[,1],0.1) | x [,1] >= quantile ( x[,1],0.9)  | x [,2] <= quantile ( x[,2],0.1) | x [,2] >= quantile ( x[,2],0.9) 
	n <- c ( sum ( ! out ), sum ( out ) )
	
	param$mean[,1] <- colMeans ( x[ ! out,] )
	param$mean[,2] <- param$mean[,1]
	param$var[,,1] <- cov ( x [ ! out,] )
	if ( uncorr == TRUE ) 
	    { 
	    param$var[1,2,1] <- 0
	    param$var[2,1,1] <- 0
	    }
	param$var[,,2] <- lambda^2 * param$var[,,1] 
	}
    
     if (! is.null(hyper_prior_mean) )
	{
	hyper_param$mean <- hyper_prior_mean 
	param$mean[,1] <- hyper_param$mean
	param$mean[,2] <- hyper_param$mean
	param$var[,,1] <- prior_scale / ( prior_df + 3 )

	if ( uncorr == TRUE ) 
	    { 
	    param$var[1,2,1] <- 0
	    param$var[2,1,1] <- 0
	    }
	param$var[,,2] <- lambda^2 * param$var[,,1] 
	
	for ( j in 1:m )
	    {
	    df <- 2 + hyper_prior_df[j]
	    sigma2 <- ( hyper_prior_df[j] * hyper_prior_scale[j] + sum ( ( param$mean[j,] - hyper_param$mean[j] )^2 ) ) / df
	    hyper_param$var[j,j] <- 1 / rgamma ( 1 , shape = df / 2  , scale = 2 / ( df * sigma2 ) )
	    } 

	#Update z
	log_p1 <- dmvnorm ( x , mean = param$mean[,2] , sigma = param$var[,,2] , log = TRUE ) 
	log_p0 <- dmvnorm ( x , mean = param$mean[,1] , sigma = param$var[,,1] , log = TRUE )  
	
	if ( q == 0 ) { theta <- rep ( 0 , n ) }
	if ( q != 0 ) { theta <- 1 / ( ( ( 1 - q ) / q ) * exp ( log_p0 - log_p1 ) + 1  )  }

	z <- rbinom ( nb_ind , size = 1 , prob = theta )
	out <- as.logical ( z )
	n <- c ( sum ( ! out ), sum ( out ) )
	}

     ## Gibbs sampling
     for ( iter in 1 : niter )
	  {

	  if ( verbose & iter %% 500 == 0 ) 
	       { 
	       writeLines ( paste ( "Iteration" , iter ) ) 
	       }

	  # Update mean and variance for normal and outlier individuals
	  # If priors are given
	  if (! is.null ( hyper_prior_mean ) )
	    {
	    for ( i in 1:2 )
		{
		if ( i == 1 ) { subset <- ! out }
		if ( i == 2 ) { subset <- out }

		mat <- diag( 1 / diag ( hyper_param$var ) )
		val <- solve ( mat + n[i] * solve ( param$var[,,i] ) )

		if ( n[i] > 1 ) { mu <- mat %*% hyper_param$mean + solve ( param$var[,,i] ) %*% colSums ( x [ subset ,] ) }
		if ( n[i] == 1 ) { mu <- mat %*% hyper_param$mean + solve ( param$var[,,i] ) %*% x [ subset ,] }
		if ( n[i] == 0 ) { mu <- mat %*% hyper_param$mean }
		
		param$mean[,i] <- rmvnorm ( 1 , mean =  val %*% mu , sigma = val )
		}

	    if ( uncorr == FALSE )
		{
		df <- n[1] + prior_df
		if ( n[1] >= 1 ) { sigma2 <- prior_scale + ( t(x[!out,])-param$mean[,1] ) %*% t ( t(x[!out,]) - param$mean[,1] ) }
		if ( n[1] == 0 ) { sigma2 <- prior_scale }
		param$var[,,1] <- riwish ( df , sigma2 )
		}

	    if ( uncorr == TRUE )
		{
		for ( j in 1:m )
		    {
		    df <- n[1] + prior_df[j]
		    if ( n[1] >= 1 ) { sigma2 <- ( prior_df[j] * prior_scale[j,j] + sum ( ( x[ !out,j] - param$mean[j,1])^2 ) ) / df }
		    if ( n[1] == 0 ) { sigma2 <- ( prior_df[j] * prior_scale[j,j] ) / df }
		    param$var[j,j,1] <- 1 / rgamma ( 1 , shape = df / 2 , scale = 2 / ( df * sigma2 ) )
		    }
		}
	    }

	  # If priors are not given (uninformative priors)
	  if ( is.null ( hyper_prior_mean ) )
	    {
	    for ( i in 1:2 )
		{
		if ( i == 1 ) { subset <- ! out }
		if ( i == 2 ) { subset <- out }

		if ( n[i] > 1 ) { mu <- colSums ( x [ subset,] ) / n[i] }
		if ( n[i] == 1 ) { mu <- x [ subset,] / n[i] }
		if ( n[i] == 0 ) { stop ( "Iteration without inlier or outlier identified: please provide some priors so that all distributions are well defined" ) }

		sigma <- param$var[,,i] / n[i]
		param$mean[,i] <- rmvnorm ( 1 , mean = mu, sigma = sigma ) 
		}

	    if ( uncorr == FALSE )
		{
		df <- n[1]
		if ( n[1] >= 1 ) { sigma2 <-  ( t(x[!out,])-param$mean[,1] ) %*% t ( t(x[!out,]) - param$mean[,1] ) }
		if ( n[1] == 0 ) { stop ( "Iteration without inlier identified: please provide some priors so that all distributions are well defined" ) }
		param$var[,,1] <- riwish ( df , sigma2 ) 
		}

	    if ( uncorr == TRUE )
		{
		for ( j in 1:m )
		    {
		    df <- n[1]
		    if ( n[1] >= 1 ) { sigma2 <- sum ( ( x[!out,j]-param$mean[j,1])^2) / df }
		    if ( n[1] == 0 ) { stop ( "Iteration without inlier identified: please provide some priors so that all distributions are well defined" ) }
		    param$var[j,j,1] <- 1 / rgamma ( 1 , shape = df/2 , scale = 2 / ( df * sigma2) )
		    }
		}
	    }
	  

	  param$var[,,2] <- lambda^2 * param$var[,,1] 


	  #Update hyper parameters if priors have been given
	  if ( ! is.null ( hyper_prior_mean ) )
	      {
	      for ( j in 1:m )
		  {
		  var <- 1 / ( 1/hyper_prior_var[j] + 2/hyper_param$var[j,j])
		  mu <- ( hyper_prior_mean[j] /hyper_prior_var[j] + sum( param$mean[j,] ) / hyper_param$var[j,j] ) * var
		  hyper_param$mean[j] <- rnorm ( 1 , mean = mu, sd = sqrt ( var ) )

		  df <- hyper_prior_df[j] + 2
		  sigma2 <- ( hyper_prior_df[j] * hyper_prior_scale[j] + sum ( ( param$mean[j,] - hyper_param$mean[j] )^2 ) ) / df
		  hyper_param$var[j,j] <- 1 / rgamma ( 1 , shape = df / 2  , scale = 2 / ( df * sigma2 ) )
		  } 
	      }

	  #Update z
	  log_p1 <- dmvnorm ( x , mean = param$mean[,2] , sigma = param$var[,,2] , log = TRUE ) 
	  log_p0 <- dmvnorm ( x , mean = param$mean[,1] , sigma = param$var[,,1] , log = TRUE )  

	  if ( q == 0 ) { theta <- rep ( 0 , n ) }
	  if ( q != 0 ) { theta <- 1 / ( ( ( 1 - q ) / q ) * exp ( log_p0 - log_p1 ) + 1  )  }

	  z <- rbinom ( nb_ind , size = 1 , prob = theta )
	  out <- as.logical(z)
	  n <- c ( sum ( !out ) , sum ( out) )

	  #Update q if prior is given
	  if ( ! is.null ( alpha ) ) { q <- rbeta ( 1 , shape1 = n[2] + alpha , shape2 = n[1] + beta ) }

	  if ( iter > burnin )
	       {
	       counts <- counts + z 
	       sum_mean <- sum_mean + param$mean[,1]
	       sum_var <- sum_var + param$var[,,1]
	       }

	  }

     post <- counts / (niter-burnin)
     outlier <- which ( post >= 0.5 )
     inlier <- which ( post < 0.5 )

     # Compute posterior mean and posterior variance for "normal" individuals
     mean_post_mean <- sum_mean / (niter-burnin)
     mean_post_var <- sum_var / (niter-burnin)

     # if verbose=TRUE, display the number of outliers identified on the screen
     if ( verbose )
	  {
	  if ( length ( outlier ) <= 1 ) { text <- "outlier" }
	  if ( length ( outlier ) > 1 ) { text <- "outliers" }
	  writeLines ( paste ( "\n" ,length ( outlier ) , text , "\n\n" ) )
	  }

     if ( length ( outlier ) != 0 ) { group [ outlier ] <- 1 }
     group [ inlier ] <- 0

     return ( list ( x = x_ini , group = group , posterior = post , lambda = lambda , post_mean = mean_post_mean , post_var = mean_post_var , standardize = standardize, inlier = inlier , outlier = outlier ) )

     }

