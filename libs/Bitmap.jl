# ********************************************************************
# SaveBitmap
#       pixMap: 256 color level image arrar
#   colorTable: 256 color scale
# ********************************************************************
function SaveBitmap( name,
  pixMap::Array{UInt8,2},
  colorTable::Array{UInt32,1} )

  # Get image dimensions
  dim = size( pixMap )
  if length( dim ) != 2
    return 0 
  end
  
  # Compute row padding
  padding = ( 4 - dim[1] & 0x3 ) & 0x3

  # Compute file dimensions. Header = 14 + 40 + ( 256 * 4 ) = 1078
  offset   = 1078
  imgBytes = dim[2] * ( dim[1] + padding )
  
  # Create file
  stream = open( name, "w" )
  
  # Save file header
  write( stream, UInt16( 0x4D42 ) )
  write( stream, UInt32[ offset + imgBytes, 0 , offset ] )

  # Save info header
  write( stream, UInt32[ 40, dim[1], dim[2], 0x80001, 0 ] )
  write( stream, UInt32[ imgBytes, 0, 0, 256, 0 ] )
  
  # Save color table
  write( stream, colorTable )
  if length( colorTable ) < 256
    fixTable = zeros( UInt32, 256 - length( colorTable ) )
    write( stream, fixTable )
  end
  
  # Save image pixels
  if padding == 0
    for i = 1:dim[2] 
      write( stream, pixMap[:,i] )
    end
  else
    zeroPad = zeros( UInt8, padding )
    for i in 1:dim[2]
      write( stream, pixMap[:,i] )
      write( stream, zeroPad )
    end
  end
  
  # Close file
  close( stream )
  
end


# ********************************************************************
# FindMass
#   massVector: mz Vector sorted in ascending order
#         mass: target mz value
#    tolerance: bi-axial tolerance to find mass 
# ********************************************************************
function FindMass( massVector, mass, tolerance )

  index  = Int( 0 )
  lower  = Int( 1 )
  higher = Int( length( massVector ) )
  
  while lower <= higher 
  
    # Compute mid element
    index = ( lower + higher ) ÷ 2

    # Go to lower portion ?
    if massVector[ index ] > ( mass + tolerance )
      higher = index - 1
      continue
    end
    
    # Go to Higher portion?
    if massVector[ index ] < ( mass - tolerance )
      lower = index + 1
      continue
    end
    
    # mass was found
    return index
  end
  
  # mass was not found
  return 0 
end


# ********************************************************************
# GetSlice
#       imzML: imzML data
#        mass: mz target value
#   tolerance: bi-axial tolerance to find mass
# ********************************************************************
function GetSlice( imzML, mass, tolerance )

  # Alloc space for slice
  width  = maximum( imzML[1,:] )
  height = maximum( imzML[2,:] )
  image  = zeros( Float64, width, height )

  for i in 1:size( imzML )[2]
    index = FindMass( imzML[3,i], mass, tolerance )
    if index != 0
     image[ imzML[1,i], imzML[2,i] ] = imzML[4,i][index]
    end
  end
  
  return image
  
end


# ********************************************************************
# IntQuant, Discretize image amplitude in 0:255 range
#   slice: Image matriz in measured arbitrary units
# ********************************************************************
function IntQuant( slice )

  # Compute scale factor for amplitude discretization
  lower  = minimum( slice )
  scale  = 255 / maximum( slice )
  dim    = size( slice )
  image  = zeros( UInt8, dim[1], dim[2] )
  
  for i in 1:length( slice )
    image[i] = convert( UInt8, floor( slice[i] * scale + 0.5 ) )
  end
  
  return image
  
end



# ********************************************************************
# TrIQ, Discretize image amplitude in 0:255 range, grouping outliers
# in the highest bin
#   slice: Image matrix in measured arbitrary units
# ********************************************************************
function GetOutlierThres( slice, prob )

  # Get bin width
  low = minimum( slice )
  upp = maximum( slice )

  # Compute histogram's bin count & separation
  if ( upp - low + 1 ) >= 100 
    bins    = 100
    step    = ( upp + (upp-low)/(bins-1) - low ) / bins
  else
    bins   = ceiling( upp ) - floor( low ) + 1
    step   = 1
  end

  # Compute histogram
  histCount = zeros( bins, 1 )
  for k in slice
    index = convert( Int64, floor( ( k - low ) / step ) + 1 )
    histCount[ index ] += 1
  end

  # Get bin index, that accounts for desired probability
  prob  = 0.95
  delta = cumsum( histCount, dims=1 ) / sum( histCount ) .- prob
  index = findmin( broadcast( abs, delta ) )[2][1]

  # Find max intensity image value 
  key   = low + step * index
  upp   = -1

  for k in slice[:]
    if key > k  &&  k > upp
      upp = k        
    end
  end
  
  return [ low, upp ]
    
end


# ********************************************************************
# Discretize image's intensity using a given intensity range 
#    slice: Image matrix in measured arbitrary units
#   bounds: max and min intensity discretizing range
#    depth: number of intensity discrete steps in final image 
# ********************************************************************
function SetPixelDepth( slice, bounds, depth )

  # Compute intensity bins
  bins   = depth - 1
  limits = collect( 
    range( bounds[1], stop = bounds[2], length = depth ) )

  # Reserve memory for output image
  imgBytes = zeros( UInt8, size( slice ) )
  nPixels  = length( slice )
  step     = limits[2] - limits[1]
  
  # Set intensity depth
  for k in 1:nPixels
    
    # Find bin
    i::Int64 = floor( ( slice[k] - limits[1] ) / step ) + 1
    if i < bins
      imgBytes[k] = i + ( slice[k] > limits[i] ) - 1
    else
      imgBytes[k] = bins
    end 

  end
  
  return imgBytes 
end


# ********************************************************************
# Discretize image's intensity removing intensity outliers 
#   slice: Image matrix in measured arbitrary units
#   depth: Number of intensity discrete steps in final image
#    prob: Proportion of pixesl to take into account in discretization
# ********************************************************************
function TrIQ( slice, depth, prob = 0.98 )
  
  return  SetPixelDepth(
    slice,
    GetOutlierThres( slice, prob ),
    depth
  )  

end

# cd( "C:/Users/LABI/Documents/_LABI/IR3/Julia" )
# img = zeros(UInt8, 32, 32)
# img[17:32,1:16] .= 1
# img[1:16,17:32] .= 1

# SaveBitmap( "test.bmp", img, UInt32[0, 0x00FFFFFF] )

