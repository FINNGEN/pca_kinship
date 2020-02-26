aberrant.plot <-
function ( x )
     {
     
     if ( missing ( x ) ) { stop ("x must be specified") }

     if ( ! is.list ( x ) || is.null ( names ( x ) ) || names ( x ) != c ( "x" , "group" , "posterior" , "lambda" , "post_mean" , "post_var" , "standardize" , "inlier" , "outlier") )
	  {
	  stop ("x must be an output from function aberrant")
	  }

     x_ini <- x$x
     n <- nrow ( x_ini )
     m <- ncol ( x_ini )

     inlier <- x$inlier
     outlier <- x$outlier
     lambda <- x$lambda
     post <- x$posterior
     post_mean <- x$post_mean
     post_var <- x$post_var
     standardize <- x$standardize
     

     if ( standardize ) { x <- scale ( x_ini ) }
     if ( ! standardize ) { x <- x_ini }
     
     if ( length ( outlier ) > 1 ) { main_title <- paste ( length ( outlier ) , "outliers" ) }
     if ( length ( outlier ) <= 1 ) { main_title <- paste ( length ( outlier ) , "outlier" ) }
     

     x_lim <- c ( min ( x[,1] ) , max ( x[,1] ) )
     y_lim <- c ( min ( x[,2] ) , max ( x[,2] ) ) 

     if ( ! is.null ( colnames ( x ) ) ) 
	    {
	    x_lab <- colnames ( x ) [1]
	    y_lab <- colnames ( x ) [2] 
	    }

     if ( is.null ( colnames ( x ) ) ) 
	    {
	    x_lab <- "1st summary statistic"
	    y_lab <- "2nd summary statistic"
	    }
      
     col_inliers <- densCols ( x [ inlier, ] ,  colramp = colorRampPalette ( c ( "grey70" , "grey60" , "grey50" , "grey40" , "grey30" , "grey20" , "grey10" , "black" ), space = "Lab") )
     col_ini <- rev ( heat.colors ( 91 ) )[ - c (1:40) ]
     col_outliers <- col_ini [ ( round ( post[outlier] , 2 ) * 100 ) - 49 ]

     if ( length ( inlier ) != 0 & length ( inlier ) != 1)
	    {
	    plot ( x [ inlier , ] , pch = 20 , xlim = x_lim , ylim = y_lim , xlab = x_lab , ylab = y_lab , main = main_title , col = col_inliers ) 
	    }
     if ( length ( inlier ) == 1)
	    {
	    plot ( x [ inlier , 1],x [ inlier , 2] , pch = 20 , xlim = x_lim , ylim = y_lim , xlab = x_lab , ylab = y_lab , main = main_title , col = col_inliers ) 
	    }
     if ( length ( outlier ) != 0 & length ( outlier ) != 1) { points ( x [ outlier , ] , col = col_outliers , pch = 20 ) }
     if ( length ( outlier ) == 1) { points ( x [ outlier,1 ] ,  x [outlier, 2 ], col = col_outliers , pch = 20 ) }
     points (  ellipse ( post_var , centre = post_mean, level = 0.99) , col = "grey30" , lty = 2 , type = "l" )
     }

