using Documenter
using julia_mzML_imzML

makedocs(
    sitename = "julia_mzML_imzML",
    format   = Documenter.HTML(),
    modules  = [julia_mzML_imzML],
		pages    = [ "Home" => "index.md" ] )


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
