# ********************************************************************
# Activate imzML package
# ********************************************************************

using Pkg

import Pkg; Pkg.add("Libz")
import Pkg; Pkg.add("Plots")

using Libz
Pkg.activate( "../../julia_mzML_imzML" )
using julia_mzML_imzML

# Data repository
#Robert Winkler. (2023). mzML mass spectrometry and imzML mass spectrometry imaging test data [Data set]. 
#Zenodo. https://doi.org/10.5281/zenodo.10084132

#The script assumes the following directory structure

#rob@uga-itx ~> tree julia* -L 1
#julia_example-data
#├── AP-SMALDI_mouse_urinary_bladder
#└── mzML
#julia_mzML_imzML
#├── docs
#├── LICENSE
#├── Manifest.toml
#├── Project.toml
#├── README.md
#├── src
#├── test
#└── timing


# ********************************************************************
# Timing mzML
# ********************************************************************

# Image load
data_dir = "../../julia_example-data/mzML"
mzML     = [ "Col_1.mzML", "Cytochrome_C.mzML", "T9_A1.mzML" ]

function MzmlTime( fileName, reps )

	spectra = []
	
	@time begin
		# Load mzML file
		for i = 1:reps
			spectra  = LoadMzml( fileName )
		end
	end
	
	return spectra
	
end


MzmlTime( joinpath( data_dir, mzML[1] ), 10 )



# ********************************************************************
# Timing imzML
# ********************************************************************

# Auxiliar function for testing image load
function SaveSlice( path, prefix, spectra, mz, tolerance, divisor )

  for k in 1:length(mz)

    # Build file name
    name = ( path * "/" * prefix * "_" *
		  string( mz[k] ? divisor, pad = 3 ) * "-" *
      string( mz[k] % divisor ) * ".bmp" )

    # Image render & save  
    SaveBitmap( name,
      IntQuant( GetSlice( spectra, mz[k] / divisor, tolerance ) ),
      ViridisPalette )
  end

end


# Image load
cd( "C:/Users/LABI/Desktop" )
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ"
imzML    = [ 
  "AP_SMALDI/HR2MSI mouse urinary bladder S096.imzML", 
	"DESI_MSI/80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML", 
	"LAESI_MSI/180817_NEG_Thaliana_Leaf_bottom_1_0841.imzML", 
	"LTP_MSI/ltpmsi-chilli.imzML" ]

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
SaveSlice( destDir, "Carcinoma", spectra, [88555], 0.005, 100 )

@time spectra = ImzmlTime( joinpath( data_dir, imzML[3] ), 10 )
SaveSlice( pww(), "Arabidopsis", spectra, [2090, 4360, 4471], 0.05, 10 )

@time spectra = ImzmlTime( joinpath( data_dir, imzML[4] ), 10 )
SaveSlice( destDir, "Chili", spectra, [621, 841, 3061], 0.1, 10 )

