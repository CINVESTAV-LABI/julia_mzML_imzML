using Libz

# Load source code
code_dir = "./"
include( joinpath( code_dir, "mzML.jl" ) )
include( joinpath( code_dir, "Decode64.jl" ) )

# Read imzML file
data_dir = "./"
scans    = LoadOffsets( joinpath( data_dir, "T9_A1.mzML" ) )

# Plot first tree scans
using Plots

plot( scans[:,1], scans[:,2] )
plot( scans[:,1], scans[:,3] )
plot( scans[:,1], scans[:,4] )
