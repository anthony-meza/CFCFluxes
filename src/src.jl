
maindir = "/Users/anthonymeza/Library/CloudStorage/OneDrive-MassachusettsInstituteofTechnology/Documents/GitHub/CFCSolubility"
datadir(x) = maindir  * "/data/" * x
plotsdir(x) = maindir  * "/plots/" * x

include("gas_flux_functions.jl")
include("example_scenarios.jl")