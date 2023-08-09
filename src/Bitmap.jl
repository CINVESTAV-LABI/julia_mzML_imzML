# ********************************************************************
# SaveBitmap
# ********************************************************************
"""
    `SaveBitmap( name, pixMap, colorTable )`

Save a discretized imzML slice as a bitmap file

# Arguments
* `      name`: Full path name of target bitmap
* `    pixMap`: [UInt8] Bidimensional matrix with image info
* `colorTable`: [UInt32] Vector with RGB colors for each gray level

# Examples
```julia
# Saves a bitmap with black and white alternate image squares
img = zeros(UInt8, 32, 32)
img[17:32,1:16] .= 1
img[1:16,17:32] .= 1
SaveBitmap( "test.bmp", img, UInt32[0, 0x00FFFFFF] )
```

The following image will be created on your hard disk

![](./Test.bmp)
"""
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
"""
    `GetSlice( imzML, mass, tolerance )`

Extract a mz-image from an imzML image data array loaded with
 `LoadImzml` function. The resulting image array preserves the same
 data type of the y-axis vectors in the file.

# Arguments
* `    imzML`: Image array loaded with `LoadImzml` function
* `     mass`: mz value for image extraction
* `tolerance`: Maximum lateral tolerance for a valid mz search

# Examples
```julia
# Load DESI MSI Carcinoma  885.55 mz slice
spectra = LoadImzml( "80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML" )
slice   = GetSlice( spectra, 885.55, 0.005 )
```
"""
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
"""
    `IntQuant( slice )`

Zero Memory intensity quantizer. Discretize the image?s continuos
 intensity range in discrete bins from 0 to 255 gray levels.
 Returns a UInt8 bidimensional matrix 

# Arguments
* `slice`: slice returned by the `GetSlice` function

# Examples
```julia
# Save DESI MSI Carcinoma 885.55 mz slice as a bitmap file
spectra = LoadImzml( "80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML" ) 
slice   = GetSlice( spectra, 885.55, 0.005 
SaveBitmap( "Slice.bmp", IntQuant( slice ), ViridisPalette )
```

The following image will be created on your hard disk

![](./Slice.bmp)
"""
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
"""
    `TrIQ( slice, depth, prob = 0.98 )`

TrIQ intensity quantizer as described in _DOI 10.7717/peerj-cs.585_.
 Discretize the image?s continuos intensity range in discrete bins
 from 0 to 255 gray levels, groupping intensity outiers in the
 highehst bin. Returns a UInt8 bidimensional matrix 

# Arguments
* `slice`: slice returned by the `GetSlice` function
* `depth`: Number ob discrete bins for intensity quantization
* ` prob`: Cummulative distribute function cutting value

# Examples
```julia
# Save DESI MSI Carcinoma 885.55 mz slice as a bitmap file
spectra = LoadImzml( "80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML" ) 
slice   = GetSlice( spectra, 885.55, 0.005 ) 
SaveBitmap( "TrIQ.bmp", TrIQ( slice, 256, 0.95 ), ViridisPalette )
```

The following image will be created on your hard disk

![](./TrIQ.bmp)
"""
function TrIQ( slice, depth, prob = 0.98 )

  return  SetPixelDepth(
    slice,
    GetOutlierThres( slice, prob ),
    depth
  )

end


# ********************************************************************
# Viridis color palette 256 color levels
# ********************************************************************
ViridisPalette = UInt32[
  0x440154, 0x440256, 0x450457, 0x450559,
  0x46075A, 0x46085C, 0x460A5D, 0x460B5E,
  0x470D60, 0x470E61, 0x471063, 0x471164,
  0x471365, 0x481467, 0x481668, 0x481769,
  0x48186A, 0x481A6C, 0x481B6D, 0x481C6E,
  0x481D6F, 0x481F70, 0x482071, 0x482173,
  0x482374, 0x482475, 0x482576, 0x482677,
  0x482878, 0x482979, 0x472A7A, 0x472C7A,
  0x472D7B, 0x472E7C, 0x472F7D, 0x46307E,
  0x46327E, 0x46337F, 0x463480, 0x453581,
  0x453781, 0x453882, 0x443983, 0x443A83,
  0x443B84, 0x433D84, 0x433E85, 0x423F85,
  0x424086, 0x424186, 0x414287, 0x414487,
  0x404588, 0x404688, 0x3F4788, 0x3F4889,
  0x3E4989, 0x3E4A89, 0x3E4C8A, 0x3D4D8A,
  0x3D4E8A, 0x3C4F8A, 0x3C508B, 0x3B518B,
  0x3B528B, 0x3A538B, 0x3A548C, 0x39558C,
  0x39568C, 0x38588C, 0x38598C, 0x375A8C,
  0x375B8D, 0x365C8D, 0x365D8D, 0x355E8D,
  0x355F8D, 0x34608D, 0x34618D, 0x33628D,
  0x33638D, 0x32648E, 0x32658E, 0x31668E,
  0x31678E, 0x31688E, 0x30698E, 0x306A8E,
  0x2F6B8E, 0x2F6C8E, 0x2E6D8E, 0x2E6E8E,
  0x2E6F8E, 0x2D708E, 0x2D718E, 0x2C718E,
  0x2C728E, 0x2C738E, 0x2B748E, 0x2B758E,
  0x2A768E, 0x2A778E, 0x2A788E, 0x29798E,
  0x297A8E, 0x297B8E, 0x287C8E, 0x287D8E,
  0x277E8E, 0x277F8E, 0x27808E, 0x26818E,
  0x26828E, 0x26828E, 0x25838E, 0x25848E,
  0x25858E, 0x24868E, 0x24878E, 0x23888E,
  0x23898E, 0x238A8D, 0x228B8D, 0x228C8D,
  0x228D8D, 0x218E8D, 0x218F8D, 0x21908D,
  0x21918C, 0x20928C, 0x20928C, 0x20938C,
  0x1F948C, 0x1F958B, 0x1F968B, 0x1F978B,
  0x1F988B, 0x1F998A, 0x1F9A8A, 0x1E9B8A,
  0x1E9C89, 0x1E9D89, 0x1F9E89, 0x1F9F88,
  0x1FA088, 0x1FA188, 0x1FA187, 0x1FA287,
  0x20A386, 0x20A486, 0x21A585, 0x21A685,
  0x22A785, 0x22A884, 0x23A983, 0x24AA83,
  0x25AB82, 0x25AC82, 0x26AD81, 0x27AD81,
  0x28AE80, 0x29AF7F, 0x2AB07F, 0x2CB17E,
  0x2DB27D, 0x2EB37C, 0x2FB47C, 0x31B57B,
  0x32B67A, 0x34B679, 0x35B779, 0x37B878,
  0x38B977, 0x3ABA76, 0x3BBB75, 0x3DBC74,
  0x3FBC73, 0x40BD72, 0x42BE71, 0x44BF70,
  0x46C06F, 0x48C16E, 0x4AC16D, 0x4CC26C,
  0x4EC36B, 0x50C46A, 0x52C569, 0x54C568,
  0x56C667, 0x58C765, 0x5AC864, 0x5CC863,
  0x5EC962, 0x60CA60, 0x63CB5F, 0x65CB5E,
  0x67CC5C, 0x69CD5B, 0x6CCD5A, 0x6ECE58,
  0x70CF57, 0x73D056, 0x75D054, 0x77D153,
  0x7AD151, 0x7CD250, 0x7FD34E, 0x81D34D,
  0x84D44B, 0x86D549, 0x89D548, 0x8BD646,
  0x8ED645, 0x90D743, 0x93D741, 0x95D840,
  0x98D83E, 0x9BD93C, 0x9DD93B, 0xA0DA39,
  0xA2DA37, 0xA5DB36, 0xA8DB34, 0xAADC32,
  0xADDC30, 0xB0DD2F, 0xB2DD2D, 0xB5DE2B,
  0xB8DE29, 0xBADE28, 0xBDDF26, 0xC0DF25,
  0xC2DF23, 0xC5E021, 0xC8E020, 0xCAE11F,
  0xCDE11D, 0xD0E11C, 0xD2E21B, 0xD5E21A,
  0xD8E219, 0xDAE319, 0xDDE318, 0xDFE318,
  0xE2E418, 0xE5E419, 0xE7E419, 0xEAE51A,
  0xECE51B, 0xEFE51C, 0xF1E51D, 0xF4E61E,
  0xF6E620, 0xF8E621, 0xFBE723, 0xFDE725 ]

