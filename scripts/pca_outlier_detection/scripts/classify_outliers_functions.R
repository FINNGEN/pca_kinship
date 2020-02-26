basic_plot_pcs <- function( data, colour_by,eigen_values = NULL,shape=NULL, n_first=4, PC_prefix="PC",...)
{
  pc_cols <- paste(PC_prefix, seq(20), sep="")   
  plots = list()
  if (! is.null(eigen_values))
  {
    var_data <- data.frame(pc=seq(length(eigen_values)), eigen_vals=eigen_values)
    var_data$var_explained <- eigen_values / sum(var_data$eigen_vals)  
    p <- ggplot(var_data) + geom_line( aes(y=var_explained, x=pc)) + scale_x_continuous( breaks=seq(1,20)) + ggtitle( "Scree plot. ")
    plot(p)
    plots[[ length(plots) + 1]] <- p
  }
  
  for( i in seq(n_first))
  {
    j <- i+1
    p <- ggplot(data ) + geom_density(aes_string(x=paste(PC_prefix,i, sep=""),colour=colour_by) )  
    plot(p)
    plots[[ length(plots) + 1]] <- p
    
    while( j <= n_first && j != i)
    {
      pc1 <- paste(PC_prefix, i, sep="")
      pc2 <- paste(PC_prefix, j, sep="")
      
      if( is.null(shape))
      {
        p <- ggplot(data) + geom_point(aes_string(x=pc1,y=pc2, colour=colour_by), size=1, ... ) + guides(colour = guide_legend(override.aes = list(size=6))) 
        plots[[ length(plots) + 1]] <- p
        plot(p)
        p <- ggplot(data) + geom_density2d(aes_string(x=pc1,y=pc2, colour=colour_by) ) 
        plot(p)
        plots[[ length(plots) + 1]] <- p
      }
      else
      {
        p <- ggplot(data) + geom_point(aes_string(x=pc1,y=pc2, colour=colour_by, shape=shape), ... ) + guides(colour = guide_legend(override.aes = list(size=6))) 
      
        plot(p)
        plots[[ length(plots) + 1]] <- p
        p <- ggplot(data) + geom_density2d(aes_string(x=pc1,y=pc2, colour=colour_by) ) 
        plot(p)
        plots[[ length(plots) + 1]] <- p
      }
      
      j <- j +1
    }
  }
  
}
