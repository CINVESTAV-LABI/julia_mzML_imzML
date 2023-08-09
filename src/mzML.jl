# *******************************************************************
# Load scans from mzML file
#   fileName: Full name path
# *******************************************************************
"""
    `LoadMzml( fileName )`

Load the mzML file as a Matrix. Each column stores a single scan,
 x-axis and y-axis scan data are stored in the first and second
 row respectively.

# Arguments
* `fileName`: Full path name of the mzML file

# Examples
```julia
spectra = LoadMzml( "T9_A1.mzML" )
size( spectra )
(2, 10)
using Plots
plot( spectra[1,4], spectra[2,4] )  
```

The following image will appear in an auxiliary window

![](./Plots.bmp)
"""
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


# *******************************************************************
# Decoding base64 blocks
# *******************************************************************

# Decoding exchange table
base64dec = Vector{UInt8}( [
  62, # (01) + -> 43d, 0x2B
   0, # (02) ,
   0, # (03) -
   0, # (04) .
  63, # (05) /
  52, # (06) 0
  53, # (07) 1
  54, # (08) 2
  55, # (09) 3
  56, # (10) 4
  57, # (11) 5
  58, # (12) 6
  59, # (13) 7
  60, # (14) 8
  61, # (15) 9
   0, # (16) :
   0, # (17) ;
   0, # (18) <
   0, # (19) =
   0, # (20) >
   0, # (21) ?
   0, # (22) @
   0, # (23) A
   1, # (24) B
   2, # (25) C
   3, # (26) D
   4, # (27) E
   5, # (28) F
   6, # (29) G
   7, # (40) H
   8, # (41) I
   9, # (42) J
  10, # (43) K
  11, # (44) L
  12, # (45) M
  13, # (46) N
  14, # (47) O
  15, # (48) P
  16, # (49) Q
  17, # (40) R
  18, # (41) S
  19, # (42) T
  20, # (43) U
  21, # (44) V
  22, # (45) W
  23, # (46) X
  24, # (47) Y
  25, # (48) Z
   0, # (49) [
   0, # (50) \
   0, # (51) ]
   0, # (52) ^
   0, # (53) _
   0, # (54) `
  26, # (55) a
  27, # (56) b
  28, # (57) c
  29, # (58) d
  30, # (59) e
  31, # (60) f
  32, # (61) g
  33, # (62) h
  34, # (63) i
  35, # (64) j
  36, # (65) k
  37, # (66) l
  38, # (67) m
  39, # (68) n
  40, # (69) o
  41, # (70) p
  42, # (71) q
  43, # (72) r
  44, # (73) s
  45, # (74) t
  46, # (75) u
  47, # (76) v
  48, # (77) w
  49, # (78) x
  50, # (79) y
  51, # (80) z
] )

# base64 decoding
function DecodeTriplet( data, base64 )

  # Character to binary equivalent
  data    = data .- 0x2A
  data[1] = base64[ data[1] ]
  data[2] = base64[ data[2] ]
  data[3] = base64[ data[3] ]
  data[4] = base64[ data[4] ]

  # Get Triplet
  decode    = Array{UInt8}( undef,3 )
  decode[1] = data[1] << 2 + data[2] >> 4
  decode[2] = data[2] << 4 + data[3] >> 2
  decode[3] = data[3] << 6 + data[4]
  return decode

end


function Decode64( data )

  # Reserve uninitialized memory
  block  = divrem( sizeof( data ), 4 )
  decode = Array{UInt8}( undef, 3*block[1] )

  # Recover binary from Chars
  index = 1
  first = 1
  last  = 3

  for i in 1:block[1]
    decode[first:last] = DecodeTriplet( data[index:index+3], base64dec )
    index += 4
    first += 3
    last  += 3
  end

  # Adjust block size
  if data[end-1] == 0x3D
    decode = decode[1:end-2]
  elseif data[end] == 0x3D
    decode = decode[1:end-1]
  end

  return decode

end
