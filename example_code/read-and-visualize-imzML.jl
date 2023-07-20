using Libz

# Load source code
code_dir = "C:/Users/LABI/Documents/_LABI/IR3/imzML_Julia/imzML"
include( joinpath( code_dir, "imzML.jl" ) )
include( joinpath( code_dir, "Decode64.jl" ) )
include( joinpath( code_dir, "Bitmap.jl" ) )
include( joinpath( code_dir, "Viridis.jl" ) )

# Load stream from file
srcPath = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/"
srcName = "AP_SMALDI/HR2MSI mouse urinary bladder S096"

# Load imzML file
hFile = open( srcPath * srcName * ".imzML" )
hIbd  = open( srcPath * srcName * ".ibd" )
wrk   = ReadImgHeader( hFile, hIbd )
close( hIbd )
close( hFile )

# Save slice image
slice = GetSlice( wrk, 798.54, 0.005 )    # Bladder
img   = IntQuant( slice )                 # Scale intensity to 256 Gray Scale
SaveBitmap( "Bladder.bmp", img, ViridisPalette )