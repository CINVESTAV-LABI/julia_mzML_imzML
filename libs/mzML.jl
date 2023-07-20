# *******************************************************************
# Load scans from mzML file
#   fileName: Full name path
# *******************************************************************
function LoadMzml( fileName )

  # Open file hFile = open( joinpath( data_dir, "Col_1.mzML" ), "r+" )
  stream    = open( fileName, "r+" )
	intensity = []
  
  # Load footer
  seekend( stream )
  skip( stream, -160 )
  raw    = read( stream, 160 )
  footer = String( raw )
  
  # Search "indexListOffset" tag
  index = match( r"<indexListOffset>(?<Ofst>\d+)</indexListOffset>", footer )
  
  if index !== nothing
  
    # Jump to index table
    seek( stream, parse( Int64, index.captures[1] ) )
  
		# Find first "index" tag
		index = FindTag( stream, r"<index +name=\"spectrum" )
		
		# Load offset list and load spectra
		if index !== nothing
			intensity = LoadSpectra( stream, ParseOffsetList( stream ) )
		end
		
	end
	
  # Close file handle
  close( stream )
  
  # Return spectral matrix
  return intensity
  
end


# *******************************************************************
# Parse indexList offsets
#   stream: Source file stream
# *******************************************************************
function ParseOffsetList( stream )

  # List of 8 elements
  entry  = zeros( Int64, 8 )
  used   = 0
  unused = 8
  
  while true
   
    # Get offset string
    raw   = readline( stream )
    index = match( r">(?<Ofst>\d+)<", String( raw ) )
    
    # Exit condition
    if index === nothing
      break
    end
    
    # Save index
    used   += 1
    unused -= 1
    entry[used] = parse( Int64, index.captures[1] )
    
    # Grow list when needed
    if unused == 0
      append!( entry, zeros( Int64, 8 ) )
      unused = 8
    end
    
  end
  
  return entry[1:used]
  
end


# *******************************************************************
# Read vector data
#   stream: Source file stream
# *******************************************************************
function ReadVector( stream, axis )

	# Read base64 vector
	base64Vec = FindTag( stream, r"<binary>([^<]+)<" )
  base64Vec = Vector{ UInt8 }( base64Vec.captures[1] )
  
  # Fill a IO stream with decoded data
	if axis.Packed == 1
		io = IOBuffer( Libz.inflate( Decode64( base64Vec ) ) )
	else
		io = IOBuffer( Decode64( base64Vec ) )
	end
	
  # Convert IO stream content as dataType vector 
  nElem = Int32( io.size / sizeof( axis.Format ) )
  out   = Array{ axis.Format }( undef, nElem )
  read!( io, out )

	return out
  
end


# *******************************************************************
# Get the decoding options for the first vector in spectrum
#    stream: Valid file stream
#   lastTag: Last tag to skip i.e. "binaryDataArray"
# *******************************************************************
function AxisConfig( stream )

  offset = position( stream )
  skip   = 0
  while true

    # Load line
    raw   = readline( stream )
    tag   = match( r"^\s*<([^\s]+)", raw )
    start = length( tag.match ) + 1

    # Is a "cvParam" tag? 
    if tag.captures[1] == "cvParam"

      # Get accesion field value
      attribute = GetAttribute( raw[ start:end ], "accession" )

      # Subtract file name length from skip field
      if attribute.captures[1] != "MS:1000796"
        continue
      else
        start += attribute.offset + length( attribute.match )
        value  = GetAttribute( raw[ start:end ], "value" )
        skip  += length( value.captures[1] )
        continue
      end

    # Exit condition 
    elseif tag.captures[1] == "binaryDataArray"
      value  = GetAttribute( raw[ start:end ], "encodedLength" )
      offset = position( stream ) - offset
      front  = ConfigureSpecDim( stream )
      front.Skip += offset - length( value.captures[1] ) - skip
      return ( front )

    end
  end
end


# *******************************************************************
# Load Spectra and return a matrix
#   stream: Source file stream
# *******************************************************************
function LoadSpectra( stream, offset )

  # Move file pointer towards first spectrum tag
  seek( stream, offset[1] )

  # Decodes axis information and skip counter
  firstAxis  = AxisConfig( stream )
  secondAxis = AxisConfig( stream )

  # Allocates memory for spectra storage
	specCount = length( offset )
	spectrum  = Array{ Any }( undef, ( 2, specCount ) )
	index     = 1
  order     = 1 + ( firstAxis.Axis == 2 )  

	for i in 1:specCount

		# Load first vector
    seek( stream, offset[ i ] )
    skip( stream, firstAxis.Skip )
		spectrum[ order, i ] = ReadVector( stream, firstAxis )

    # Load second vector
    skip( stream, secondAxis.Skip )
		spectrum[ xor( order,3 ), i ] = ReadVector( stream, secondAxis )
	end

  return( spectrum )
  
end


