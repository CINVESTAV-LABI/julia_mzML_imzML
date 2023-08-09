module julia_mzML_imzML

# external dependencies
using Libz

# mzML Functions
export LoadMzml

# imzML Functions
export LoadImzml

# Bitmap functions
export GetSlice
export IntQuant
export TrIQ
export SaveBitmap

# Viridis color scheme
export ViridisPalette

# Source files
include( "Common.jl" )
include( "mzML.jl" )
include( "imzML.jl" )
include( "Bitmap.jl" )

end
