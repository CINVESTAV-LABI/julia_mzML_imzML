module LabiMzLib

# external dependencies
using Libz

# mzML Functions
export LoadMzml

# imzML Functions
export LoadImzml

# bitmap Functions
export GetSlice
export IntQuant
export SaveBitmap

# viridis color scheme
export ViridisPalette


include( "Bitmap.jl" )
include( "Decode64.jl" )
include( "Viridis.jl" )
include( "mzML.jl" )
include( "imzML.jl" )

end
