using Libz
code_dir = "."
include( joinpath( code_dir, "Decode64.jl" ) )
include( joinpath( code_dir, "Common.jl" ) )
include( joinpath( code_dir, "mzML.jl" ) )

# Read imzML file
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/imzML_Julia/mzML"
spectra  = LoadMzml( joinpath( data_dir, "T9_A1.mzML" ) )

# LC-MS
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/mzML/LC-MS/Col"
spectra  = LoadMzml( joinpath( data_dir, "Col_1.mzML" ) )

# EsiProt
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/mzML/ESIProt"
spectra  = LoadMzml( joinpath( data_dir, "Cytochrome_C.mzML" ) )

using Plots
plot( spectra[1,1], spectra[2,1] )
plot( spectra[1,2], spectra[2,2] )
plot( spectra[1,3], spectra[2,3] )
plot( spectra[1,4], spectra[2,4] )
plot( spectra[1,32], spectra[2,32] )




# ********************************************************************
# imzML TEST
# ********************************************************************

# Load source code
srcPath  = "C:/Users/LABI/.julia/dev/LabiMzLib/src"
include( joinpath( srcPath, "Decode64.jl" ) )
include( joinpath( srcPath, "Common.jl" ) )
include( joinpath( srcPath, "imzML.jl" ) )

include( joinpath( srcPath, "Viridis.jl" ) )
include( joinpath( srcPath, "Bitmap.jl" ) )

# Auxiliar function for testing image load
function SaveSlice( path, prefix, spectra, mz, tolerance, divisor )

  for k in 1:length(mz)

    # Build file name
    name = ( path * "/" * prefix * "_"
        * string( mz[k] ÷ divisor, pad = 3 ) * "-" 
        * string( mz[k] % divisor ) * ".bmp" )

    # Image render & save  
    SaveBitmap( name,
      IntQuant( GetSlice( spectra, mz[k] / divisor, tolerance ) ),
      ViridisPalette )
  end

end

# Timing
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/AP_SMALDI"
fileName = "HR2MSI mouse urinary bladder S096.imzML"
@time LoadSpectraImg( joinpath( data_dir, fileName ) )


# Path for image storage
destDir  = "C:/Users/LABI/Desktop/Test"

# AP_SMALDI: Mouse Bladder
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/AP_SMALDI"
fileName = "HR2MSI mouse urinary bladder S096.imzML"
spectra  = LoadSpectraImg( joinpath( data_dir, fileName ) )
SaveSlice( destDir, "Mouse", spectra, [74153, 74354, 79854], 0.005, 100 )

# LAESI_MSI: Arabidopsis
# https://doi.org/10.5281/zenodo.3678472
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/LAESI_MSI"
fileName = "180817_NEG_Thaliana_Leaf_bottom_1_0841.imzML"
spectra  = LoadSpectraImg( joinpath( data_dir, fileName ) )
SaveSlice( destDir, "Arabidopsis", spectra, [2090, 4360, 4471], 0.05, 10 )

# https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100131/ColAd_Individual/ColAd_Individual.zip
# DESI_MSI: Carcinoma
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/DESI_MSI"
fileName = "40TopL,10TopR,30BottomL,20BottomR-centroid.imzML"
spectra  = LoadSpectraImg( joinpath( data_dir, fileName ) )
SaveSlice( destDir, "Carcinoma40", spectra, [88555], 0.005, 100 )

data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/DESI_MSI"
fileName = "80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML"
spectra  = LoadSpectraImg( joinpath( data_dir, fileName ) )
SaveSlice( destDir, "Carcinoma80", spectra, [88555], 0.005, 100 )

data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/DESI_MSI"
fileName = "120TopL, 90TopR, 110BottomL, 100BottomR-centroid.imzML"
spectra  = LoadSpectraImg( joinpath( data_dir, fileName ) )
SaveSlice( destDir, "Carcinoma120", spectra, [88555], 0.005, 100 )


# LTP-MSI: Chili
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/LTP_MSI"
fileName = "ltpmsi-chilli.imzML"
spectra  = LoadSpectraImg( joinpath( data_dir, fileName ) )
SaveSlice( destDir, "Chili", spectra, [621, 841, 3061], 0.1, 10 )



# ********************************************************************
# TrIQ test
# ********************************************************************
data_dir = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/DESI_MSI"
fileName = "80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML"
spectra  = LoadSpectraImg( joinpath( data_dir, fileName ) )

slice    = GetSlice( spectra, 885.55, 0.005 )
destDir  = "C:/Users/LABI/Desktop/Test"
SaveBitmap( joinpath( destDir, "TrIQ.bmp" ),
  TrIQ( slice, 256, 0.95 ),
  ViridisPalette )

# ********************************************************************
# DOWNLOAD SAMPLE FILES
# ********************************************************************

samplesDir = "C:/Users/LABI/Desktop/Download"

# AP_SMALDI: Mouse Bladder
download(
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2014/11/PXD001283/HR2MSImouseurinarybladderS096.imzML",
  joinpath( samplesDir, "HR2MSImouseurinarybladderS096.imzML" ) )

download(
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2014/11/PXD001283/HR2MSImouseurinarybladderS096.ibd",
  joinpath( samplesDir, "HR2MSImouseurinarybladderS096.ibd" ) )
  
# LAESI_MSI: Arabidopsis
download(
  "https://zenodo.org/record/3678473/files/180817_NEG_Thaliana_Leaf_bottom_1_0841.imzML",
  joinpath( samplesDir, "180817_NEG_Thaliana_Leaf_bottom_1_0841.imzML.imzML" ) )

download(
  "https://zenodo.org/record/3678473/files/180817_NEG_Thaliana_Leaf_bottom_1_0841.ibd",
  joinpath( samplesDir, "180817_NEG_Thaliana_Leaf_bottom_1_0841.ibd" ) )

# DESI_MSI: Carcinoma
download(
  "https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100131/ColAd_Individual/ColAd_Individual.zip",
  joinpath( samplesDir, "ColAd_Individual.zip" ) )  
 


