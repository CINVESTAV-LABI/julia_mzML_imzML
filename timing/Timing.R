library( MALDIquant )
library( MALDIquantForeign )


LoadingTime <- function ( fileName ) {

  elapsed <- vector( 'integer', 10 )

  for( i in 1:10 ) {
    start_time <- Sys.time()
    spectra <- MALDIquantForeign::importImzMl( fileName, removeEmptySpectra = FALSE )
    elapsed[i] <- Sys.time() - start_time
  }

  return( elapsed )
  
}

# Data repository
#Robert Winkler. (2023). mzML mass spectrometry and imzML mass spectrometry imaging test data [Data set]. 
#Zenodo. https://doi.org/10.5281/zenodo.10084132

data_dir <- "../../julia_example-data/imzML_AP_SMALDI"
fileName <- "HR2MSImouseurinarybladderS096.imzML"
time     <- LoadingTime( file.path( data_dir, fileName ) )

# Ubuntu 23.04 16GB, Intel® Celeron® N5105 × 4 [1] 1.284737 1.241130 1.212906 1.251302 1.285321 1.288718 1.329186 1.317943 1.452790 1.275222
