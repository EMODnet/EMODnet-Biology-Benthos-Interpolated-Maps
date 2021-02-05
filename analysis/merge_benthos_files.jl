include("../scripts/BenthosInterp.jl")
using NCDatasets

# Set directories:
# * where the results are stored
resultdir = "../product/netCDF/1-UniformL/"
# * where the output will be written
outputdir = "../product/netCDF/1-UniformL/Combined/"
outputfile = "Benthos_combined_V1.nc"

# Set domain and resolution
# (could be read from one of the result files instead of hardcoding)
domain = [-16., 9., 45., 66.]; # [West East South North]
Δlon = 0.1
Δlat = 0.1
longrid = domain[1]:Δlon:domain[2]
latgrid = domain[3]:Δlat:domain[4]

# Generate list of files to be merged
# (by default: only netCDF files)
resultfilelist = filter(x->endswith(x, ".nc"), readdir(resultdir));
nfiles = length(resultfilelist)
@info("Working on $(nfiles) files")

# Create a new (empty) netCDF file that will stored all the species
isdir(outputdir) ? @debug("Directory already exists") : mkdir(outputdir)
@info("Creating new netCDF file")
create_nc_results_merged(joinpath(outputdir, outputfile), longrid, latgrid)

# Open the merged file for editing and loop on the individual files
NCDatasets.Dataset(joinpath(outputdir, outputfile), "a") do nc
    for (ispec, resfile) in enumerate(resultfilelist)
        @info("Reading results from file $(resfile)")

        lon, lat, field, error, domain, scientificname, aphiaID = read_results(joinpath(resultdir, resfile))
        @info(size(field))
        nc["probability"][:,:,ispec] = field
        nc["probability_error"][:,:,ispec] = error
        nc["aphiaID"][ispec] = aphiaID
        nc["speciesname"][ispec] = scientificname
    end
end
