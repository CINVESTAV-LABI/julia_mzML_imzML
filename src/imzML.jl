# include( "Common.jl" );
# *******************************************************************
# Load Spectra and return a matrix
#   fileName: Full name path
# *******************************************************************
"""
    LoadImzml( fileName )

Load an imzML file as a matrix. Each column stores x-pixel position,
 y-pixel position, x-axis data and y-axis data.

# Arguments
* `fileName`: Full path name of the imzML file

# Examples
```julia
# Load DESI MSI Carcinoma image data
spectra = LoadImzml( "80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML" )
size( spectra )
(4, 18632)
```
"""
function LoadImzml( fileName )

  # Open file handles
  fileName = split( fileName, "." )[1]
  stream   = open( fileName * ".imzML" )
  hIbd     = open( fileName * ".ibd" )

  # Get axes types and image dimensions
  axis   = AxesConfigImg( stream )
  imgDim = GetImgDimensions( stream )
  format = [ axis[1].Format, axis[2].Format ]

  # Locate spectrum attributes
  start  = position( stream )  
  attr   = GetSpectrumAttributes( stream, hIbd ) 

  # Load spectra
  seek( stream, start )
  spectra = LoadImgData( stream, hIbd, attr, imgDim[3], format )

  close( stream )
  close( hIbd )
  return spectra

end


# *******************************************************************
# Get axes value type
# *******************************************************************
function AxesConfigImg( stream )

  # Locate which axis is defined at first
  tag   = FindTag( stream, r"^\s*<(referenceableParamGroup )" )
  value = GetAttribute( tag.captures[1], "intensityArray" )
  order = 1 + ( value !== nothing )

  # Read first axis configuration
  axis  = Array{ SpecDim, 1 }( undef, 2 )
  axis[ order ] = ConfigureSpecDim( stream )

  # Read second axis configuration
  FindTag( stream, r"^\s*<(referenceableParamGroup )" )
  axis[ xor( order,3 ) ] = ConfigureSpecDim( stream )
  return axis

end


# *******************************************************************
# Get vector's storage options and image dimensions
# *******************************************************************
function GetImgDimensions( stream )

  # Looks for "scanSettings" tag
  # FindImgTag( stream, "scanSettings" )
  FindTag( stream, r"^\s*<(scanSettings )" )

  # Initial values for dimension retrieve
  n         = 2
  dim       = [ 0, 0, 0 ]
  currLine  = ""
  matchInfo = RegexMatch
  
  while( n > 0 )
  
    # Next field
    currLine  = readline( stream )
    matchInfo = match( r"^\s*<(cvParam)", currLine )
    index     = length( matchInfo.captures[1] )
    matchInfo = GetAttribute( currLine[index:end], "accession" )
    
    # Get axis identity
    if( matchInfo.captures[1] == "IMS:1000042"        # max X
    ||  matchInfo.captures[1] == "IMS:1000043" )      # max Y
    
      # Read dimension's pixels 
      axis      = matchInfo.captures[1][end] - '1'
      index    += matchInfo.offsets[1] + length( matchInfo.captures[1] )
      matchInfo = GetAttribute( currLine[index:end], "value" )
      dim[axis] = parse( Int32, matchInfo.captures[1] )
      n        -= 1
    end
  end

  # Load stored spectra counter
  matchInfo = FindTag( stream, r"^\s*<spectrumList(.+)" )
  matchInfo = GetAttribute( matchInfo.captures[1], "count" )
  dim[ 3 ]  = parse( Int32, matchInfo.captures[1] )
  
  return dim
end


# *******************************************************************
# Length of "spectrum" tag without attribute values
# *******************************************************************
function  GetSpectrumTag( stream )

  offset = position( stream )
  tag    = FindTag( stream, r"^\s*<spectrum (.+)" )
  first  = 1

  while true
    
    value = match( r"[^=]+=\"([^\"]+)\"", tag.captures[1][first:end] )
    if value === nothing
      break
    end
    first  += value.offsets[1] + length( value.captures[1] )
    offset += length( value.captures[1] )
  end

  return offset
end

# seek( stream, 7141 )

# *******************************************************************
# Locate spectrum attribute's 
#   [1]  1: x position stored first   2: y position stored first
#   [2]  Second stored position
#   [3]  1: mzArray stored first      2: intensityArray stored first
#   [4]  Second stored axis
#   [5]  skip chars to find first dimension
#   [6]  skip chars to find second dimension
#   [7]  skip chars to find vector length
#   [8]  skip chars to find next spectrum  
# *******************************************************************
function GetSpectrumAttributes( stream, hIbd )

  # Reserve memory for vector configuration
  skip   = Vector{ UInt32 }( undef, 8 )
  offset = GetSpectrumTag( stream )

  # Get axis order and position of pixel coordinate values
  tag     = FindTag( stream, r" accession=\"IMS:100005(\d)\"(.+)" )
  skip[1] = tag.captures[1][1] + 1 - '0'
  skip[2] = xor( skip[1], 3 )

  value   = match( r" value=\"(\d+)\".+", tag.captures[2] )
  skip[5] = position( stream ) - offset - length( value.match ) - 2
  value   = match( r" value=\"\d+\".+", readline( stream ) )
  skip[6] = value.offset

  # Get axis order 
  offset  = position( stream )
  tag     = FindTag( stream, r"^\s*<referenceableParamGroupRef(.+)" )
  value   = GetAttribute( tag.captures[1], "ref" )
  skip[3] = ( value.captures[1] == "intensityArray" ) + 3
  skip[4] = xor( skip[3], 7 )

  k = 2
  while k != 0

    tag = FindTag( stream, r" accession=\"IMS:100010(\d)\"(.+)" )

    # Set IBD offset for first spectrum
    if tag.captures[1][1] == '2'
      value = GetAttribute( tag.match, "value" )
      seek( hIbd, parse( Int32, value.captures[1] ) )
      k -= 1
      continue

    elseif tag.captures[1][1] == '3'
      skip[7]  = position( stream ) - offset - length( tag.match )
      offset   = position( stream )
      k -= 1
    end
  end

  # Compute characters to skip 
  FindTag( stream, r"^\s*</spectrum" )
  skip[8] = position( stream ) - offset
  return skip

end


# *******************************************************************
# Read spectral data
# *******************************************************************
function LoadImgData( stream, hIbd, attr, imgDim, format )

  # Reserve spectra memory & update IBD seek offset
  spectra = Array{Any}( undef, ( 4, imgDim ) )
  contador = 0;

  for k in 1:imgDim

    # Load image coordinates
    skip( stream, attr[5] )
    value                = FindTag( stream, r"value=\"(\d+)\"" )
    spectra[ attr[1],k ] = parse( Int32, value.captures[1] )

    skip( stream, attr[6] )
    value                = FindTag( stream, r"value=\"(\d+)\"" )
    spectra[ attr[2],k ] = parse( Int32, value.captures[1] )

    # Get vector length
    skip( stream, attr[7] )
    value   = FindTag( stream, r"value=\"(\d+)\"" )
    nPoints = parse( Int32, value.captures[1] )

    # Reserve memory and read spectrum
    spectra[ 3,k ] = Array{ format[1] }( undef, nPoints )
    spectra[ 4,k ] = Array{ format[2] }( undef, nPoints )
    read!( hIbd, spectra[ attr[3],k ] )
    read!( hIbd, spectra[ attr[4],k ] )
    skip( stream, attr[8] )

  end

  return spectra

end
