# ********************************************************************
# mzML TEST
# ********************************************************************

# Load libraries
using Libz
using Plots

# Path for mzML libraries
code_dir = "../libs"

# Path for image storage
destDir  = "~/Test"

include( joinpath( code_dir, "Decode64.jl" ) )
include( joinpath( code_dir, "Common.jl" ) )
include( joinpath( code_dir, "mzML.jl" ) )

# Single mzML file
data_dir = "/home/rob/julia_mzML_imzML_testdata/mzML_single"
spectra_1  = LoadMzml( joinpath( data_dir, "T9_A1.mzML" ) )
 
# LC-MS
data_dir = "/home/rob/julia_mzML_imzML_testdata/mzML_LC-MS"
spectra_2 = LoadMzml( joinpath( data_dir, "Col_1.mzML" ) )

plot( spectra_2[1,32], spectra_2[2,32] )

# ESIprot
data_dir = "/home/rob/julia_mzML_imzML_testdata/mzML_ESIprot"
spectra_3 = LoadMzml( joinpath( data_dir, "Cytochrome_C.mzML" ) )

#plot( spectra_3[1,32], spectra_3[2,32] )


