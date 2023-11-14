library( viridis )
library( MALDIquant )
library( MALDIquantForeign )

#Change to R script directory
#setwd( "~/julia_mzML-imzML/timing" )

#TrIQ algorithm, see
#Rosas-Román I, Winkler R. Contrast optimization of mass spectrometry imaging (MSI) data visualization by threshold intensity quantization (TrIQ).
#PeerJ Comput Sci. 2021 Jun 9;7:e585. doi: 10.7717/peerj-cs.585. PMID: 34179452; PMCID: PMC8205298.
source( "TrIQ.R" )

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

ImzmlLoadingTime <- function ( fileName, reps ) {

  elapsed <- vector( 'integer', reps )

  for( i in 1:reps ) {
    start_time <- Sys.time()
    spectra <- MALDIquantForeign::importImzMl( fileName, removeEmptySpectra = FALSE )
    elapsed[i] <- Sys.time() - start_time
  }
	
	# Display results
	print( "The loading time is:" )
	print( elapsed )
	
	# Save TrIQ
	
  return( spectra )
  
}

data_dir <- "../../julia_example-data"
imzML    <- c( 
  "imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML",
  "imzML_DESI/ColAd_Individual/80TopL, 50TopR, 70BottomL, 60BottomR/80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML",
  "imzML_LA-ESI/180817_NEG_Thaliana_Leaf_bottom_1_0841.imzML",
  "imzML_LTP/ltpmsi-chilli.imzML")

# SMALDI Mouse urinary bladder
specList <- ImzmlLoadingTime( file.path( data_dir, imzML[1] ), 10 )
imgList  <- GetSlice( c( 741.53, 743.54, 798.54 ), specList, 0.005 )
dim( imgList ) <- c(3, 1)
PlotSlices( GlobalTrIQ( imgList, 256, 0.98 ), 256 )

# DESI Carcinoma
specList <- ImzmlLoadingTime( file.path( data_dir, imzML[2] ), 10 )
imgList  <- GetSlice( c( 885.55 ), specList, 0.005 )
dim( imgList ) <- c(1, 1)
PlotSlices( GlobalTrIQ( imgList, 256, 0.98 ), 256 )

# LAESI Arabidopsis
specList <- ImzmlLoadingTime( file.path( data_dir, imzML[3] ), 10 )
imgList  <- GetSlice( c( 209, 436, 447.1 ), specList, 0.05 )
dim( imgList ) <- c(3, 1)
PlotSlices( GlobalTrIQ( imgList, 64, 0.98 ), 64 )

# LTP Chili
specList <- ImzmlLoadingTime( file.path( data_dir, imzML[4] ), 10 )
imgList  <- GetSlice( c( 62.1, 84.1, 306.1 ), specList, 0.05 )
dim( imgList ) <- c(3, 1)
PlotSlices( GlobalTrIQ( imgList, 128, 0.98 ), 128 )
