# julia_mzML_imzML

## Local Installation

1. Download the package from <https://github.com/CINVESTAV-LABI/julia_mzML_imzML>.

2. Download and Install the following libraries: `Libz`, `Plots` package running the following script. this code only needs to be run the first time you run this example 

    ```julia
    import Pkg; Pkg.add("Libz")
    import Pkg; Pkg.add("Plots")
   ```

3. Activate `julia_mzML_imzML` package running the following script

   

   ```julia
   using Libz
   using Pkg
   Pkg.activate( "C:/your/dowwnload/folder/julia_mzML_imzML" )
   using julia_mzML_imzML
   ```

   `C:/your/dowwnload/folder` is the path on your computer, where the Github repository was downloaded. Now your are able to execute the test scripts.



## Loading mzML files

1. Make sure you have a mzML file. In this example, we provide a script for downloading available public mzML files.

   ``````julia
   samplesDir = "C:/some/data/folder/"
   
   # LC-ESI MS: Arabidopsis
   download(
     "https://zenodo.org/record/8185092/files/Col_1.mzML?download=1",
     joinpath( samplesDir, "Col_1.mzML" ) ) 
   
   # Replace Col_1 with Cytochrome_C to download ESI-MS Cytochrome C file
   # Replace Col_1 with T9_A1 to downloadw LTP-MS Arabidopsis file
   ``````

2. Load your data in Julia

   ```julia
   # Load mzML file
   spectra  = LoadMzml( joinpath( samplesDir, "Col_1.mzML" ) )
   ```

   Now the scans are loaded in the vector `spectra`. Each row corresponds to a single scan, the first column contains a vector with the x-axis and its corresponding y-axis is stored in the second column.

3. Plot a scan. In the following example we plot the fourth scan of `Col_1.mzML` file

   ```julia
   # Plot scan
   using Plots
   plot( spectra[1,4], spectra[2,4] )
   ```
   ![](.\test\mzML.png)



## Loading imzML files

1. The following example loads the *DESI_MSI Carcinoma*  imzML file from a public repository.

   ```julia
   samplesDir = "C:/some/data/folder/"
   
   # DESI_MSI: Carcinoma
   download(
     "https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100131/ColAd_Individual/ColAd_Individual.zip",
     joinpath( samplesDir, "ColAd_Individual.zip" ) )  
   ```

2. Extract the **imzML & ibd files** `80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML` is located inside the folder `80TopL, 50TopR, 70BottomL, 60BottomR` **Note: There are many other samples and subfolders in the dataset** 

3. Load the imzML file in memory

   ```julia
   # AP_SMALDI: Mouse Bladder
   fileName = "80TopL, 50TopR, 70BottomL, 60BottomR-centroid.imzML"
   spectra  = LoadImzml( joinpath( samplesDir, fileName ) )
   ```

4. Extract a mz slice

   ```julia
   # Extract image slice 
   slice = GetSlice( spectra, 885.55, 0.005 )
   ```

5. Render the image as bitmap

   ```julia
   # Save slice as bitmap with Zero-Memory quantizier
   SaveBitmap( joinpath( samplesDir, "Slice.bmp" ),
     IntQuant( slice ),
     ViridisPalette )
   ```

   The following image will be created in your `samplesDir` folder

   ![](.\test\Slice.bmp)

6. You can improve the dynamic range of the image using the TrIQ algorithm.

   
   ```julia
   # Improve dynamic range with TrIQ alorithm
   SaveBitmap( joinpath( samplesDir, "TrIQ.bmp" ),
     TrIQ( slice, 256, 0.95 ),
     ViridisPalette )  
   ```

   ![](.\test\TrIQ.bmp)

 # Example data

Robert Winkler. (2023). mzML mass spectrometry and imzML mass spectrometry imaging test data [Data set]. 
Zenodo. <https://doi.org/10.5281/zenodo.10084132>
