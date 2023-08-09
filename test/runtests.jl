using LabiMzLib
using Plots
using Test

@testset "LabiMzLib.jl" begin

  # mzML Test
	data_dir = "C:/Users/LABI/Documents/_LABI/IR3/imzML_Julia/mzML"
	scans    = LoadMzml( joinpath( data_dir, "T9_A1.mzML" ) )
		
	# Plot first scan
	plot( scans[:,1], scans[:,2] )
	
	# Load stream from file
	srcPath = "C:/Users/LABI/Documents/_LABI/IR3/(Data)/TrIQ/"
	srcName = "AP_SMALDI/HR2MSI mouse urinary bladder S096"

	# Load imzML file
	hFile = open( srcPath * srcName * ".imzML" )
	hIbd  = open( srcPath * srcName * ".ibd" )
	wrk   = LoadImzml( hFile, hIbd )
	close( hIbd )
	close( hFile )


end
