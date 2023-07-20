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

data_dir <- "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/AP_SMALDI"
fileName <- "HR2MSI mouse urinary bladder S096.imzML"
time     <- LoadingTime( file.path( data_dir, fileName ) )

   R         Julia 
 80.83069  12.67549
 81.38165
 89.50832
115.58657
122.19773
 99.82199
 83.32184
 85.32897
 87.17396
 97.87446