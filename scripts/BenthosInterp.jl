usecartopy = true

using GridInterpolations
using Interpolations
using Missings

"""
    read_coords_species(datafile, species)

Read the coordinates for presence and abscence events stored in the file
`datafile` for the selected species.

## Examples
```julia-repl
julia> lon_pre, lat_pre, lon_abs, lat_abs = read_data_phyto("specs4Diva.csv")
```
"""
function read_coords_species(datafile::String, species::Union{String,SubString{String}}, occurtype=1)
    data = readdlm(datafile, ',');
    colnames = data[1,:]
    dates = data[2:end,2]
    lon = data[2:end,3]
    lat = data[2:end,4]
    col_index = findall(colnames .== species)[1]
    @info("Column index for $(species): $(col_index)");
    occur_species = data[2:end,col_index]

	if occurtype == 1
		occur_pre = findall(occur_species .== "TRUE")
    	occur_abs = findall(occur_species .== "FALSE")
	elseif occurtype == 2
		occur_pre = findall(occur_species .== "1")
    	occur_abs = findall(occur_species .== "0")
	end

    lon_presence = lon[occur_pre]
    lat_presence = lat[occur_pre]
    lon_absence = lon[occur_abs]
    lat_absence = lat[occur_abs];
    return Float64.(lon_presence), Float64.(lat_presence), Float64.(lon_absence), Float64.(lat_absence)
end

"""
    get_species_slug(filename)

Return the a species 'slug' based on the species name

## Examples
```julia-repl
julia> slug = get_species_slug("Schizomavella_(Schizomavella)_auriculata")
"Schizomavella_Schizomavella_auriculata"
```
"""
function get_species_slug(specname::String)::String
    speciesslug = replace(replace(specname, " " => "_"), "." => "")
    speciesslug = replace(speciesslug, "(" => "")
    speciesslug = replace(speciesslug, ")" => "")
    return speciesslug
end

"""
    get_species_list(datafile)

Return the list of species names stored in file `datafile`

## Examples
```julia-repl
julia> specnames = get_species_list("specs4Diva.csv")
```
"""
function get_species_list(datafile::String)::Array
    specnames = split(readline(datafile), ",")[5:end]
    return specnames
end

"""
    read_specnames(datafile)

Return the a dictionary where keys are the and the values the

## Examples
```julia-repl
julia> spec_dict = read_specnames("spe.csv")
```
"""
function read_specnames(datafile::String)::Dict
    d = Dict([])
    f = readline(datafile)
    for lines in readlines(datafile)
        d[split(lines, ",")[1]] = split(lines, ",")[2]
    end
    return d
end

"""
    plot_heatmap(longrid, latgrid, dens, lonpre, latpre, lonabs, latabs,
                 titletext, figname, vmin, vmax, usecartopy)

Create a plot of the heatmap `dens` and overlay the positions of the
presence (`lonpre`, `latpre`) and absence (`lonabs`, `latabs`) data.

If `usecartopy` is set as true, then the
[`cartopy`](https://scitools.org.uk/cartopy/docs/latest/) module is used
to create the map.

If `figname` is defined, then the plot is saved as a figure.

## Example
```julia-repl
julia> plot_heatmap(longrid, latgrid, dens, lonpre, latpre, lonabs, latabs,
             "Diplocirrus stopbowitzi", "Diplocirrus_stopbowitzi_dens.png",
             0., 1., false)
```
"""
function plot_heatmap(longrid::StepRangeLen, latgrid::StepRangeLen,
    dens::Array, lonpre::Vector, latpre::Vector, lonabs::Vector, latabs::Vector,
    titletext::String; figname::String="",
    vmin::Float64=0., vmax::Float64=1.0, usecartopy=false)

    llon, llat = ndgrid(longrid, latgrid)
    myproj = ccrs.PlateCarree()
    fig = PyPlot.figure(figsize=(12,8))
    if usecartopy
        ax = PyPlot.subplot(111, projection=myproj)
    else
        ax = PyPlot.subplot(111)
    end

    ax.plot(lonpre, latpre, "wo", markersize=1., zorder=3, alpha=.25)
    ax.plot(lonabs, latabs, "ko", markersize=1., zorder=3, alpha=.25)
    pcm = ax.pcolor(llon, llat, dens, cmap=PyPlot.cm.hot_r, zorder=2,
                    vmin=vmin, vmax=vmax)
    colorbar(pcm, orientation="vertical")
    if usecartopy
        decorate_map_domain(ax)
    end

    title(titletext)
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
        PyPlot.close()
    else
        PyPlot.show()
    end
end


"""
    plot_error(longrid, latgrid, dens, titletext, figname, usecartopy)

Create a plot of the error field

If `usecartopy` is set as true, then the
[`cartopy`](https://scitools.org.uk/cartopy/docs/latest/) module is used
to create the map.

If `figname` is defined, then the plot is saved as a figure.

## Example
```julia-repl
julia> plot_error(longrid, latgrid, cmpe,
                  "Diplocirrus stopbowitzi",
                  "Diplocirrus_stopbowitzi_error.png",
                  false)
```
"""
function plot_error(longrid::StepRangeLen, latgrid::StepRangeLen,
    error::Array, titletext::String=""; figname::String="", usecartopy=false)

    llon, llat = ndgrid(longrid, latgrid)
    fig = PyPlot.figure(figsize=(12,8))
    if usecartopy
        ax = PyPlot.subplot(111, projection=myproj)
    else
        ax = PyPlot.subplot(111)
    end

    pcm = ax.pcolor(llon, llat, error, zorder=2, cmap=PyPlot.cm.RdYlGn_r)
    colorbar(pcm, orientation="vertical")
    if usecartopy
        decorate_map_domain(ax)
    end

    title(titletext)
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
        PyPlot.close()
    else
        PyPlot.show()
    end
end

"""
    make_plot_presence_absence(lon_pre, lat_pre, lon_abs, lat_abs, species, figname)

Plot the locations of absence and presence in 2 subplots and save it as `figname`.

## Example
```julia-repl
julia> make_plot_presence_absence(lon_pre, lat_pre, lon_abs, lat_abs,
"Abra_alba_data", "data_locations.jpg")
```
"""
function make_plot_presence_absence(lon_pre::Vector, lat_pre::Vector,
	                                lon_abs::Vector, lat_abs::Vector,
									species::String=""; figname::String="",
									domain=[-16., 9., 45., 66.], dlon=4., dlat=4.,
									usecartopy=false)

    fig = PyPlot.figure(figsize=(8, 8))
	if usecartopy
		ax = PyPlot.subplot(121, projection=myproj)
	else
		ax = PyPlot.subplot(121)
	end

    ax.plot(lon_pre, lat_pre, "ko", markersize=.2)
    title("Presence of $(species)")

	if usecartopy
		decorate_map_domain(ax, true; domain=domain, dlon=dlon, dlat=dlat)
	end

	if usecartopy
		ax = PyPlot.subplot(122, projection=myproj)
	else
		ax = PyPlot.subplot(122)
	end

    ax.plot(lon_abs, lat_abs, "ko", markersize=.2)
    title("Absence")

	if usecartopy
		decorate_map_domain(ax, true; domain=domain, dlon=dlon, dlat=dlat)
	end

    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
    end

    PyPlot.close()
end

if usecartopy
    using PyCall
    ccrs = pyimport("cartopy.crs")
    gridliner = pyimport("cartopy.mpl.gridliner")
    cfeature = pyimport("cartopy.feature")
    mticker = pyimport("matplotlib.ticker")
    myproj = ccrs.PlateCarree()
    coast = cfeature.GSHHSFeature(scale="high");
    mpl = pyimport("matplotlib");
    cartopyticker = pyimport("cartopy.mpl.ticker")
    lon_formatter = cartopyticker.LongitudeFormatter()
    lat_formatter = cartopyticker.LatitudeFormatter()

    function decorate_map_domain(ax, plotcoast=true;
        domain=[-16., 9., 45., 66.], dlon::Float64=4., dlat::Float64=4.)

        PyPlot.grid(linewidth=0.2, zorder=6)
        if plotcoast
            ax.add_feature(coast, facecolor="#363636",
                edgecolor="k", zorder=5)
        end
        ax.set_xlim(domain[1], domain[2])
        ax.set_ylim(domain[3], domain[4])
        ax.set_xticks(domain[1]:dlon:domain[2])
        ax.set_yticks(domain[3]:dlat:domain[4])
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
    end
end


"""
    create_nc_results(filename, lons, lats, field, speciesname)

Write a netCDF file `filename` with the coordinates `lons`, `lats` and the
heatmap `field`. `speciesname` is used for the 'title' attribute of the netCDF.

## Examples
```julia-repl
julia> create_nc_results("Bacteriastrum_interp.nc", lons, lats, field,
    "Bacteriastrum")
```
"""
function create_nc_results(filename::String, lons, lats, field,
                           speciesname::String="";
                           valex=-999.9,
                           varname = "heatmap",
                           long_name = "Heatmap",
						   domain = [-180., 180., -90., 90.]
                           )
    Dataset(filename, "c") do ds

        # Dimensions
        ds.dim["lon"] = length(lons)
        ds.dim["lat"] = length(lats)
        #ds.dim["time"] = Inf # unlimited dimension

        # Declare variables
		nccrs = defVar(ds, "crs", Int64, ())
    	nccrs.attrib["grid_mapping_name"] = "latitude_longitude"
    	nccrs.attrib["semi_major_axis"] = 6371000.0 ;
    	nccrs.attrib["inverse_flattening"] = 0 ;

        ncfield = defVar(ds, varname, Float64, ("lon", "lat"))
        ncfield.attrib["missing_value"] = Float64(valex)
        ncfield.attrib["_FillValue"] = Float64(valex)
		ncfield.attrib["units"] = "1"
        ncfield.attrib["long_name"] = long_name
		ncfield.attrib["coordinates"] = "lat lon"
		ncfield.attrib["grid_mapping"] = "crs" ;

		# No time variable needed here
        """
        nctime = defVar(ds,"time", Float32, ("time",))
        nctime.attrib["missing_value"] = Float32(valex)
        nctime.attrib["units"] = "seconds since 1981-01-01 00:00:00"
        nctime.attrib["long_name"] = "time"
        """

        nclon = defVar(ds,"lon", Float32, ("lon",))
        # nclon.attrib["missing_value"] = Float32(valex)
        nclon.attrib["_FillValue"] = Float32(valex)
        nclon.attrib["units"] = "degrees_east"
        nclon.attrib["long_name"] = "Longitude"
		nclon.attrib["standard_name"] = "longitude"
		nclon.attrib["axis"] = "X"
		nclon.attrib["reference_datum"] = "geographical coordinates, WGS84 projection"
		nclon.attrib["valid_min"] = -180.0 ;
		nclon.attrib["valid_max"] = 180.0 ;

        nclat = defVar(ds,"lat", Float32, ("lat",))
        # nclat.attrib["missing_value"] = Float32(valex)
        nclat.attrib["_FillValue"] = Float32(valex)
        nclat.attrib["units"] = "degrees_north"
		nclat.attrib["long_name"] = "Latitude"
		nclat.attrib["standard_name"] = "latitude"
		nclat.attrib["axis"] = "Y"
		nclat.attrib["reference_datum"] = "geographical coordinates, WGS84 projection"
		nclat.attrib["valid_min"] = -90.0 ;
		nclat.attrib["valid_max"] = 90.0 ;

        # Global attributes
		ds.attrib["title"] = "$(long_name) based on presence/absence of $(speciesname)"
		ds.attrib["institution"] = "GHER - University of Liege, Deltares, VLIZ"
		ds.attrib["source"] = "spatial interpolation of presence/absence data"
        ds.attrib["project"] = "EMODnet Biology Phase III"
        ds.attrib["comment"] = "Original data prepared by Deltares"
        ds.attrib["data_authors"] = "Peter Herman (Peter.Herman@deltares.nl)"
        ds.attrib["processing_authors"] = "C. Troupin (ctroupin@uliege), A. Barth (a.barth@uliege.be)"
		ds.attrib["publisher_name"] = "VLIZ"
		ds.attrib["publisher_url"] = "http://www.vliz.be/"
		ds.attrib["publisher_email"] = "info@vliz.be"
        ds.attrib["created"] = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
		ds.attrib["geospatial_lat_min"] = domain[3]
		ds.attrib["geospatial_lat_max"] = domain[4]
		ds.attrib["geospatial_lon_min"] = domain[1]
		ds.attrib["geospatial_lon_max"] = domain[2]
		ds.attrib["geospatial_lat_units"] = "degrees_north"
		ds.attrib["geospatial_lon_units"] = "degrees_east"
		ds.attrib["license"] = "GNU General Public License v2.0"
		ds.attrib["citation"] = "A. Barth, P. Hermann and C. Troupin (2020). Probability maps for different benthos species in the North Sea."
		ds.attrib["acknowledgement"] = "European Marine Observation Data Network (EMODnet) Biology project (EASME/EMFF/2017/1.3.1.2/02/SI2.789013), funded by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund"
		ds.attrib["tool"] = "DIVAnd"
		ds.attrib["tool_version"] = "v2.6.5"
		ds.attrib["tool_doi"] = "10.5281/zenodo.4306494"
		ds.attrib["language"] = "Julia 1.5.3"
		ds.attrib["Conventions"] = "CF-1.7"
		ds.attrib["netcdf_version"] = "4"

        # Define variables
        ncfield[:] = field

        nclon[:] = lons
        nclat[:] = lats;

    end
end;

"""
    write_nc_error(filename, error)

Write the error field in the file `filename` created by `create_nc_results`

## Examples
```julia-repl
julia> write_nc_error(filename, error)
```
"""
function write_nc_error(filename::String, error::Array  ;
    valex=-999.9, varname = "heatmap",
    long_name = "Heatmap")

    Dataset(filename, "a") do ds
        ncerror = defVar(ds, varname * "_error", Float64, ("lon", "lat"))
        ncerror.attrib["missing_value"] = Float64(valex)
        ncerror.attrib["_FillValue"] = Float64(valex)
		ncerror.attrib["units"] = "1"
        ncerror.attrib["long_name"] = long_name * "_error"
		ncerror.attrib["coordinates"] = "lat lon"
		ncerror.attrib["grid_mapping"] = "crs" ;
        ncerror[:] = error
    end
end

"""
	read_substrate(datafile)

Read the coordinates and the correlation length (2D field) from the file `datafile`.

## Example
```julia-repl
julia> lon, lat, g = read_substrate("substrate_gini_impurity.nc")
```
"""
function read_substrate(datafile::String)
    NCDatasets.Dataset(datafile, "r") do nc
        lonL = nc["lon"][:]
        latL = nc["lat"][:]
        gimpurity = nc["substrate_gini_impurity"][:]

        return lonL::Array{Float64,1}, latL::Array{Float64,1}, gimpurity::Array{Float64,2}
    end
end


"""
	interp_horiz(londata, latdata, data, longrid, latgrid)

Perform a bilinear interpolation of a 2D field defined by the coordinates
`(londata, latdata)` and the values `data`, onto the grid defined by the
vectors `longrid` and `latgrid`.
The interpolation is only performed over the area where data are available,
i.e., no extrapolation is performed.

## Example
```julia-repl
julia> fieldinterp = interp_horiz(londata, latdata, data, longrid, latgrid)
```
"""
function interp_horiz(londata, latdata, data, longrid, latgrid)

    # Find the coordinates where the interpolation can be performed
    # (no extrapolation)
    goodlon = (longrid .<= londata[end]) .& (longrid .>= londata[1]);
    goodlat = (latgrid .<= latdata[end]) .& (latgrid .>= latdata[1]);

    # Create the interpolator
    itp = Interpolations.interpolate((londata, latdata), data, Gridded(Linear()))
    #fieldinterpolated = itp(longrid, latgrid);

    # Perform the interpolation
    # (only within the domain of interest)
    lon_interp = longrid[goodlon]
    lat_interp = latgrid[goodlat]
    field_interpolated = itp(lon_interp, lat_interp);

    return lon_interp, lat_interp, field_interpolated, findall(goodlon), findall(goodlat)
end


"""
    create_nc_results_merged(filename, lons, lats)

Create a netCDF file `filename` with the coordinates `lons`, `lats` and the
heatmap `field`.

## Examples
```julia-repl
julia> create_nc_results("Bacteriastrum_interp.nc", lons, lats, field,
    "Bacteriastrum")
```
"""
function create_nc_results_merged(filename::String, lons, lats, field,
                           valex=-999.9,
                           varname = "heatmap",
                           long_name = "Heatmap",
						   domain = [-180., 180., -90., 90.]
                           )
    Dataset(filename, "c") do ds

        # Dimensions
        ds.dim["lon"] = length(lons)
        ds.dim["lat"] = length(lats)
        #ds.dim["time"] = Inf # unlimited dimension

        # Declare variables
		nccrs = defVar(ds, "crs", Int64, ())
    	nccrs.attrib["grid_mapping_name"] = "latitude_longitude"
    	nccrs.attrib["semi_major_axis"] = 6371000.0 ;
    	nccrs.attrib["inverse_flattening"] = 0 ;

        ncfield = defVar(ds, varname, Float64, ("lon", "lat"))
        ncfield.attrib["missing_value"] = Float64(valex)
        ncfield.attrib["_FillValue"] = Float64(valex)
		ncfield.attrib["units"] = "1"
        ncfield.attrib["long_name"] = long_name
		ncfield.attrib["coordinates"] = "lat lon"
		ncfield.attrib["grid_mapping"] = "crs" ;

		# No time variable needed here
        """
        nctime = defVar(ds,"time", Float32, ("time",))
        nctime.attrib["missing_value"] = Float32(valex)
        nctime.attrib["units"] = "seconds since 1981-01-01 00:00:00"
        nctime.attrib["long_name"] = "time"
        """

        nclon = defVar(ds,"lon", Float32, ("lon",))
        # nclon.attrib["missing_value"] = Float32(valex)
        nclon.attrib["_FillValue"] = Float32(valex)
        nclon.attrib["units"] = "degrees_east"
        nclon.attrib["long_name"] = "Longitude"
		nclon.attrib["standard_name"] = "longitude"
		nclon.attrib["axis"] = "X"
		nclon.attrib["reference_datum"] = "geographical coordinates, WGS84 projection"
		nclon.attrib["valid_min"] = -180.0 ;
		nclon.attrib["valid_max"] = 180.0 ;

        nclat = defVar(ds,"lat", Float32, ("lat",))
        # nclat.attrib["missing_value"] = Float32(valex)
        nclat.attrib["_FillValue"] = Float32(valex)
        nclat.attrib["units"] = "degrees_north"
		nclat.attrib["long_name"] = "Latitude"
		nclat.attrib["standard_name"] = "latitude"
		nclat.attrib["axis"] = "Y"
		nclat.attrib["reference_datum"] = "geographical coordinates, WGS84 projection"
		nclat.attrib["valid_min"] = -90.0 ;
		nclat.attrib["valid_max"] = 90.0 ;

        # Global attributes
		ds.attrib["title"] = "$(long_name) based on presence/absence of $(speciesname)"
		ds.attrib["institution"] = "GHER - University of Liege, Deltares, VLIZ"
		ds.attrib["source"] = "spatial interpolation of presence/absence data"
        ds.attrib["project"] = "EMODnet Biology Phase III"
        ds.attrib["comment"] = "Original data prepared by Deltares"
        ds.attrib["data_authors"] = "Peter Herman (Peter.Herman@deltares.nl)"
        ds.attrib["processing_authors"] = "C. Troupin (ctroupin@uliege), A. Barth (a.barth@uliege.be)"
		ds.attrib["publisher_name"] = "VLIZ"
		ds.attrib["publisher_url"] = "http://www.vliz.be/"
		ds.attrib["publisher_email"] = "info@vliz.be"
        ds.attrib["created"] = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
		ds.attrib["geospatial_lat_min"] = domain[3]
		ds.attrib["geospatial_lat_max"] = domain[4]
		ds.attrib["geospatial_lon_min"] = domain[1]
		ds.attrib["geospatial_lon_max"] = domain[2]
		ds.attrib["geospatial_lat_units"] = "degrees_north"
		ds.attrib["geospatial_lon_units"] = "degrees_east"
		ds.attrib["license"] = "GNU General Public License v2.0"
		ds.attrib["citation"] = "A. Barth, P. Hermann and C. Troupin (2020). Probability maps for different benthos species in the North Sea."
		ds.attrib["acknowledgement"] = "European Marine Observation Data Network (EMODnet) Biology project (EASME/EMFF/2017/1.3.1.2/02/SI2.789013), funded by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund"
		ds.attrib["tool"] = "DIVAnd"
		ds.attrib["tool_version"] = "v2.6.5"
		ds.attrib["tool_doi"] = "10.5281/zenodo.4306494"
		ds.attrib["language"] = "Julia 1.5.3"
		ds.attrib["Conventions"] = "CF-1.7"
		ds.attrib["netcdf_version"] = "4"

        # Define variables
        ncfield[:] = field

        nclon[:] = lons
        nclat[:] = lats;

    end
end;
