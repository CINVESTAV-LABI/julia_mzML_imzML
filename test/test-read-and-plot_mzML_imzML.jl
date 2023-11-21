# ********************************************************************
# Activate imzML package
# ********************************************************************

using Pkg

import Pkg; Pkg.add("Libz")
import Pkg; Pkg.add("Plots")

using Libz
using Plots

Pkg.activate("../julia_mzML_imzML")
using julia_mzML_imzML

# Data repository
#Robert Winkler. (2023). mzML mass spectrometry and imzML mass spectrometry imaging test data [Data set]. 
#Zenodo. https://doi.org/10.5281/zenodo.10084132

#The script assumes the following directory structure
#rob@uga-itx ~> tree julia* -L 1
#julia_example-data
#├── imzML_AP_SMALDI
#├── imzML_DESI
#├── imzML_LA-ESI
#├── imzML_LTP
#└── mzML
#julia_mzML_imzML
#├── docs
#├── src
#├── test
#└── timing


# Output directory
cd( "../Downloads/" )

# ********************************************************************
# Reading and printing mzML spectra (single, multile scans, LC-MS
# ********************************************************************

# Load data

data_dir = "../julia_example-data/mzML"
mzML     = [ "Col_1.mzML", "Cytochrome_C.mzML", "T9_A1.mzML" ]

# Define output file, load spectra and plot

Test1= string(mzML[1],".pdf")
spectra  = LoadMzml( joinpath(data_dir, mzML[1] ) )
mz =  plot( spectra[1,124], spectra[2,124], seriestype=:sticks, lc=:black, legend=false )
xlabel!("m/z")
ylabel!("intensity")
hline!([0], lc=:black) 
savefig(mz, Test1)

Test1BPC= string(mzML[1],"_BPC.pdf")
mz =  plot( maximum(spectra[2,:]), lc=:blue, legend=false )
xlabel!("scan")
ylabel!("base peak intensity")
savefig(mz, Test1BPC)

Test2= string(mzML[2],".pdf")
spectra  = LoadMzml( joinpath(data_dir, mzML[2] ) )
mz =  plot( spectra[1,8], spectra[2,8], seriestype=:sticks, lc=:black, legend=false )
xlabel!("m/z")
ylabel!("intensity")
hline!([0], lc=:black) 
savefig(mz, Test2)

Test3= string(mzML[3],".pdf")
spectra  = LoadMzml( joinpath(data_dir, mzML[3] ) )
mz =  plot( spectra[1,8], spectra[2,8], seriestype=:sticks, lc=:black, legend=false )
xlabel!("m/z")
ylabel!("intensity")
hline!([0], lc=:black) 
savefig(mz, Test3)

# Help function for Image loading and saving
function SaveSlice( path, prefix, spectra, mz, tolerance, divisor )

  for k in 1:length(mz)

    # Build file name
    name = ( path * "/" * prefix * "_" *
		  string( mz[k] % divisor, pad = 3 ) * "-" *
      string( mz[k] % divisor ) * ".bmp" )

    # Image render & save  
    SaveBitmap( name,
#     IntQuant( GetSlice( spectra, mz[k] / divisor, tolerance)),
     TrIQ( GetSlice( spectra, mz[k] / divisor, tolerance ) , 256, 0.95),
      ViridisPalette )
  end

end

# ********************************************************************
# Reading and printing imzML images
# ********************************************************************

data_dir = "/home/rob/julia_example-data/"
imzML    = [ 
  "imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML", 
  "imzML_DESI/ColAd_Individual/80TopL, 50TopR, 70BottomL, 60BottomR/80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML", 
  "imzML_LA-ESI/180817_NEG_Thaliana_Leaf_bottom_1_0841.imzML", 
  "imzML_LTP/ltpmsi-chilli.imzML" ]

function ImzmlTime( fileName, reps )

	# Load mzML file
	spectra = [];
	for i = 1:reps
		spectra  = LoadImzml( fileName )
	end
	
	return spectra

end

@time spectra = ImzmlTime( joinpath( data_dir, imzML[1] ), 10 )
SaveSlice( pwd(), "Mouse", spectra, [74153, 74354, 79854], 0.005, 100 )

@time spectra = ImzmlTime( joinpath( data_dir, imzML[2] ), 10 )
SaveSlice( pwd(), "Carcinoma", spectra, [88555], 0.005, 100 )

@time spectra = ImzmlTime( joinpath( data_dir, imzML[3] ), 10 )
SaveSlice( pwd(), "Arabidopsis", spectra, [2090, 4360, 4471], 0.05, 10 )

@time spectra = ImzmlTime( joinpath( data_dir, imzML[4] ), 10 )
SaveSlice( pwd(), "Chili", spectra, [621, 841, 3061], 0.1, 10 )
