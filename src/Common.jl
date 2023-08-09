# *******************************************************************
# VectorConfig, structure with axis decoding instructions
# *******************************************************************
mutable struct SpecDim
  Format::DataType     # Float32, Float64
  Packed::Bool         # 0: none     1: zlib compression
  Axis::Int32          # 1: m/z      2: Amplitude
  Skip::Int64          # Bytes to skip
end

# *******************************************************************
# Fill VectorConfig fields from "cvParam" tags
# *******************************************************************
function ConfigureSpecDim( stream )

  # Vector configuration
  axis = SpecDim( Float64, false, 1, 0 )

  # Initial values
  offset    = position( stream )
  currLine  = ""
  matchInfo = RegexMatch
  
  while true
  
    # Next field
    currLine  = readline( stream )
    matchInfo = match( r"^\s*<(cvParam)", currLine )
    
    if matchInfo === nothing
      matchInfo = match( r"^\s*", currLine )
      axis.Skip = position( stream ) - offset - 
                  length( currLine ) + length( matchInfo.match )
      return axis
    end

    index     = length( matchInfo.captures[1] )
    matchInfo = GetAttribute( currLine[index:end], "accession" )
    
    if matchInfo.captures[1] == "MS:1000515"        # intensity array
      axis.Axis = 2
      continue
      
    elseif matchInfo.captures[1] == "MS:1000519"    # 32-bit integer
      axis.Format = Int32
      continue

    elseif matchInfo.captures[1] == "MS:1000521"    # 32-bit float
      axis.Format = Float32
      continue
      
    elseif matchInfo.captures[1] == "MS:1000522"    # 64-bit integer
      axis.Format = Int64
      continue

    elseif matchInfo.captures[1] == "MS:1000574"		# zlib compresion
			axis.Packed = true
			continue

    end
  end
end


# *******************************************************************
# Read lines in file until regex match
#   stream: Source file stream
#    regex: target tag e.g. "binaryDataArray"
# *******************************************************************
function FindTag( stream, regex )

  while true
    isTag = match( regex, readline( stream ) )
    if isTag !== nothing
      return isTag
    end
  end  
end


# *******************************************************************
# Retrieve attribute/value pair
#   source: String with atr="value" pairs
#      tag: Optional attribute, e.g. "\\spectrum"
# *******************************************************************
function GetAttribute( source, tag = "([^=]+)" )

  # Build regex and retrieve attribute/value pair
  regStr = Regex( "\\s" * tag * "=\"([^\"]*)\"" )
  return( match( regStr, source ) )
  
end


