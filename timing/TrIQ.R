# **********************************************************************
# Loading external libraries
# **********************************************************************
if( is.loaded( "viridis" ) == FALSE ) 
  library( viridis )

if( is.loaded( "MALDIquantForeign" ) == FALSE )
  library( MALDIquantForeign )


# **********************************************************************
# TrIQ Functions
# **********************************************************************

# Find intensity threshold closer to cumulative probability
#   pixMap: Numeric matrix with image intensity
#     prob: Probability, 0.95-0.99 recommended
GetOutlierThres <- function( pixMap, prob ) {
  
  # Compute bins
  lim <- range( pixMap )
  
  if( lim[2] - lim[1] + 1 >= 100 ) {
    bins <- 100
    brk  <- seq( 
      from       = lim[1], 
      to         = lim[2] + ( lim[2] - lim[1] ) / ( bins - 1 ), 
      length.out = bins+1 )
  }
  
  else {
    brk <- floor( lim[1] ):( ceiling( lim[2] ) + 1 )
  }

  # Compute histogram
  h <- hist( pixMap, 
    breaks         = brk, 
    right          = FALSE, 
    include.lowest = FALSE, 
    plot           = FALSE )

  # Get the bin closer to "prob"
  top   <- prob - cumsum( h$counts ) / sum( h$counts )
  delta <- abs( top )
  index <- 1 + which( min(delta) == delta )[1]
  
  # return range
  return( c( lim[1] , max( pixMap[ pixMap < h$breaks[index] ] ) ) )
}


# Adjust intensity dinamic range
#   pixMap: Input image
#   bounds: Intensity range on final image
#    depth: Number of grey levels
SetPixelDepth <- function( pixMap, bounds, depth ) {
  
  # Reserve memory for output matrix
  size     <- dim( pixMap )
  digital  <- matrix( 0, nrow=size[1], ncol=size[2] )
  
  # Compute bin Width
  binWidth <- ( depth - 1 ) / ( bounds[2] - bounds[1]  )
  value    <- c( 0, 0, depth - 1 )
  
  # Amplitude quantization
  for( i  in 1:length( pixMap ) ) {
    numb       <- ( pixMap[i] - bounds[1] ) * binWidth
    whole      <- trunc( numb )
    value[2]   <- whole - ( numb == whole )
    digital[i] <- value[ ( numb > 1 ) + ( pixMap[i] > bounds[2] ) + 1 ]
  }
  
  return( digital )
}


# Remove pixel intensity outliers
#   pixMap: Input image
#    depth: Number of grey levels
#     prob: Target cumulative probability
TrIQ <- function( pixMap, depth, prob=0.98 ) {

  # Compute new dinamic range
  bounds <- GetOutlierThres( pixMap, prob )
  
  # Set intensity dinamyc range
  return( SetPixelDepth( pixMap, bounds, depth ) )
}


# Apply TrIQ to a group of slices
#   slices: Slices list array
#    depth: Grey levels
#     prob: Target cumulative probability
GlobalTrIQ <- function( slices, depth, prob=0.98 ) {

  # Get intensity range
  bounds  <- range( sapply( slices, GetOutlierThres, prob ) )
  normImg <- lapply( slices, SetPixelDepth, bounds, depth )
  dim( normImg ) <- dim( slices )
  
  # Display bounds
  print( bounds )
  
  # return normalization range
  return( normImg )
}



# ****************************************************************************
# Auxiliar Functions
# ****************************************************************************

# Get slice from imzML object, set NAs to minimum value
#    mass: m/z value
#   imzML: List of MALDIquant spectrum objects
#     tol: peak tolerance
GetSlice <- function( mass, imzML, tol ) {

  # Set NA's to minimum intensity value
  FixNA <- function( mzImg ) {
 
    index <- is.na( mzImg )
    if( length( index ) ) {
      mzImg[ index ] <- min( mzImg, na.rm=TRUE )
    }
    
    return( mzImg )
  }
    
  # Load Slice

  # Array to list
  slice  <- lapply( 
    1:length( mass ), 
    function( x,y ) FixNA( 
      msiSlices( imzML, center=mass[x], tolerance=tol )[,,1] ), 
    sliceArr ) 
 
  return( slice )
}



PlotSlices <- function( slices, depth ) {
  
  # Set drawing options
  old <- par( no.readonly=TRUE )
  par( mfcol = dim( slices ), mar = rep( 0.2, 4 ) )
  
  # Draw slices
  lapply( slices, function(x) 
    image( x, asp=ncol(x)/nrow(x), axes=FALSE, col=viridis( depth ) ) )
  
  # Restore options
  par( old )
}


# Median Filter with 3x3 mask
#   pixMap: input image
MedianFilter <- function( pixMap ) {
  
  width  <- ncol(pixMap)
  height <- nrow(pixMap)
  
  target <- matrix( 0, nrow=height, ncol=width )
  
  for( j in 2:(width-1) ) {
    for( i in 2:(height-1) ) {
      index <- (i-1):(i+1)
      target[i,j] <- median( 
        c( pixMap[ index, j-1 ], 
           pixMap[ index, j ],
           pixMap[ index, j+1 ] ) )
    }
  }
  
  return( target )
}


# *************************************************************************
# Plot color scale
#  bounds: Intensity range
#   depth: Color depth
#    tick: Index for labeling
# *************************************************************************
PlotColorScale <- function( bounds, depth, tick=0 ) { 
  
  # Get label format
  exponent <- floor( log10( bounds[2] ) ) %/% 3
  if( exponent < 0 ) {
    bounds <- bounds * ( 10 ^ ( 3 * exponent ) )
  } else {
    bounds <- bounds / ( 10 ^ ( 3 * exponent ) )
  }
  
  decimals <- 1 + floor( log10( bounds[2] ) )
  region   <- ( bounds[2] - bounds[1] ) / depth
  if( ( region - floor( region ) == 0 ) ) {
    decimals <- 4
  }
    
  formatStr <- c( "%3.2f", "%3.1f", "%3.1f", "%d" )
  
  # Get color list
  scheme <- viridis( depth )
  
  # Modify margins
#  old    <- par( no.readonly=TRUE )
#  par( mar = c( 1.5,15.5,2.5,15.5) )

  # Plot color scale
  plot( c(0, 1), c(0,depth), axes=FALSE, 
    type = "n", xlab = "", xaxt="n", xlim=c(0,1),
    main = "" , ylab = "", yaxt ="n" )
  
  invisible( lapply( 0:(depth-1), function(x) 
    rect( 0, x, 100, x+1, col=scheme[ x+1 ], border=scheme[ x+1 ] ) ) )
  
  # Compute vertical axis
  if( length(tick) == 1 ) {
    tick <- GetColorDivisions( depth )
  }
  
  axis( 
    side   = 4,
    at     = tick,
    las    = 2,
    labels = sprintf( 
      formatStr[ decimals ], 
      tick * region + bounds[1] ) )
  
  # Print scale factor
  if( exponent != 0 ) {
    prefix <- parse( text=
      paste( 'x10', as.character(exponent*3), sep="^" ) )
    mtext( prefix, 3, 0.25, cex=0.8 )
  }
 
  # Restore options
#  par( old )
}


GetColorDivisions <- function( depth ) {
    
  divisor <- c(3:7)
  remain  <- depth %% divisor
  
  # Is depth a multiple of divisor ?
  if( sum( remain == 0 ) ) {
    return( seq( 
      from = 0, 
      to   = depth, 
      by   = depth / divisor[ which( remain == 0 )[1] ] ) )
  } 
    
  # Returns an approximation
  quot    <- depth %/% divisor
  index   <- which.min( remain / divisor )[1]
  return( seq( from=remain[index], to=depth, by=quot[index] )  )
}
