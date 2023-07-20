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


