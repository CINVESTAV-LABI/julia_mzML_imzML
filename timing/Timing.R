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

#Römpp A, Guenther S, Schober Y, Schulz O, Takats Z, Kummer W, Spengler B.
#ProteomeXchange dataset PXD001283. 2014.
#http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001283
#ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/11/PXD001283
#https://www.ebi.ac.uk/pride/archive/projects/PXD001283
#Publication
#DOI: 10.1002/anie.200905559, PubMed: 20397170,
#Des: Römpp A, Guenther S, Schober Y, Schulz O, Takats Z, Kummer W, Spengler B;
#Histology by mass spectrometry: label-free tissue characterization obtained
#from high-accuracy bioanalytical imaging., Angew Chem Int Ed Engl, 2010 May 17, 49, 22, 3834-8

data_dir <- "../julia_example-data/bladder"
fileName <- "HR2MSImouseurinarybladderS096.imzML"
time     <- LoadingTime( file.path( data_dir, fileName ) )
