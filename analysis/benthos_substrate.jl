using CSV
using DataFrames
using PyPlot
using Glob
using NCDatasets
using Interpolations
#using HypothesisTests
using Dates
using OceanPlot

struct Format2020Benthos
    filename::String
    df
end

function Format2020Benthos(filename)
    df = DataFrame(CSV.read(datafile; missingstring="NA",
              dateformat="yyyy-mm-ddTHH:MM:SSZ",
              truestrings=["TRUE"], falsestrings=["FALSE"]))
    return Format2020Benthos(filename,df)
end

function listnames(d::Format2020Benthos)
    return String.(propertynames(d.df)[5:end])
end

function loadbyname(d::Format2020Benthos,years,scientificname)
    df = d.df
    obstime = df.eventDate
    ids = df.eventNummer # sic
    lon = df.decimalLongitude
    lat = df.decimalLatitude

    value = getproperty(df,Symbol(scientificname))
    sel = (years[1] .<= Dates.year.(obstime)  .<= years[end]) .& (.!ismissing.(value))
    return lon[sel],lat[sel],obstime[sel],Float64.(value[sel]),ids[sel]

end


#cd(expanduser("~/src/EMODnet-Biology-Benthos-Interpolated-Maps/data"))

datadir = "../data/"
datafile = joinpath(datadir, "raw_data/specs4Diva.csv")



d = Format2020Benthos(datafile)
years = [-1000,3000]

scientificnames = listnames(d)
@show listnames(d)

longrid = -16.:0.1:9.
latgrid = 45.:0.1:66.


#scatter(lon,lat,10,Int.(value))
#colorbar()


function loadsubstrate(datadir)
    substrate_filenames = glob("substrate*nc",datadir)

    # exclude generic Seabed class
    substrate_filenames = filter(f -> basename(f) != "substrate_Seabed.nc",substrate_filenames)

    substrate_names = replace.(replace.(basename.(substrate_filenames),".nc" => ""),"substrate_" => "")
    ds = Dataset(substrate_filenames[1])
    s_lon = ds["longitude"][:]
    s_lat = ds["latitude"][:]
    close(ds)

    s_data = zeros(length(s_lon),length(s_lat),length(substrate_names))

    for i = 1:length(substrate_names)
        varname = substrate_names[i]
        ds = Dataset(substrate_filenames[i])
        s_data[:,:,i] = nomissing(ds[varname][:,:],0)
        close(ds)
    end

    s_lat = reverse(s_lat)
    s_data = reverse(s_data,dims=2)

    return s_lon,s_lat,substrate_names,s_data
end

function count_per_substrate(d,years,refvalue,s_lon,s_lat,substrate_names,s_data)
    count_per_substrate = zeros(length(scientificnames),length(substrate_names))
    fraction_per_substrate = zeros(length(scientificnames),length(substrate_names))

    for l = 1:length(scientificnames)
        scientificname = scientificnames[l]
        lon,lat,obstime,value,ids = loadbyname(d,years,scientificname)

        presence = value .== refvalue

        for i = 1:length(substrate_names)
            itp = interpolate((s_lon,s_lat), s_data[:,:,i], Gridded(Linear()))

            s_interp = itp.(lon,lat) .> 0.5
            count_per_substrate[l,i] = sum(s_interp .& presence)
        end
    end

    fraction_per_substrate = count_per_substrate ./ sum(count_per_substrate,dims=2)

    return count_per_substrate,fraction_per_substrate
end

function plot_per_substrate(scientificnames,substrate_names,frac_presence_per_substrate,figtitle)
    figure(figsize=(7,10))
    im = imshow(100*frac_presence_per_substrate)
    ax = gca()
    #ax.margins(0.5)
    ax.set_yticks(0:length(scientificnames)-1)
    ax.set_yticklabels(replace.(scientificnames,"_" => " "))
    ax.set_xticks(0:length(substrate_names)-1)
ax.set_xticklabels(replace.(substrate_names,"." => " "),rotation=90)
    clb = colorbar(im,fraction=0.046, pad=0.04)
    #clb = colorbar()
    clb.ax.set_title("%")
    title(figtitle)
    gcf().tight_layout()
end

s_lon,s_lat,substrate_names,s_data = loadsubstrate(datadir)
presence_per_substrate,frac_presence_per_substrate = count_per_substrate(d,years,1,s_lon,s_lat,substrate_names,s_data)
absence_per_substrate, frac_absence_per_substrate = count_per_substrate(d,years,0,s_lon,s_lat,substrate_names,s_data)
@show presence_per_substrate
@show absence_per_substrate

#=

mkpath("../figures")
plot_per_substrate(scientificnames,substrate_names,frac_presence_per_substrate,"Presence")
savefig("../figures/Persence.png")

plot_per_substrate(scientificnames,substrate_names,frac_absence_per_substrate,"Absence")
savefig("../figures/Absence.png")

close("all")
=#

"""
     make sure a pixel is only in one class
"""
function only_one_substrate!(s_data)
    tmp = s_data[:,:,1] .== 1

    for i = 2:size(s_data,3)
        view(s_data,:,:,i)[tmp] .= 0
        tmp = tmp .| (s_data[:,:,i] .== 1)
    end
end

only_one_substrate!(s_data)
findall( (sum(s_data,dims=3)) .> 1)

nsubstrate = size(s_data,3)
substrate_index = sum(s_data .* reshape(1:nsubstrate,(1,1,nsubstrate)),dims=3)[:,:,1]

#pcolormesh(s_data[:,:,1]')

substrate_index[substrate_index .== 0] .= NaN

#=
pcolormesh(s_lon,s_lat,substrate_index')
sb = colorbar()
colorbar()
=#
@show substrate_names


function giniimpurity_ref(y_subset,classes)
    gini_index = 1.
    probability = zeros(length(classes))

    for i = 1:length(classes)
       probability[i] = sum(y_subset .== classes[i])/length(y_subset)
    end
    gini = 1 - sum(probability.^2)
    return gini
end

function giniimpurity(y_subset,classes)
    count2 = 0
    nsubset = length(y_subset)
    @inbounds for i = 1:length(classes)
        count = 0
        for j = 1:nsubset
            if y_subset[j] == classes[i]
                count += 1
            end
        end

        count2 += count^2
    end
    gini = 1 - count2/(nsubset^2)
    return gini
end

subset = fill(1,1000)
classes = 1:10

using BenchmarkTools
@btime giniimpurity_ref(subset,classes)
@btime giniimpurity(subset,classes)
#@code_warntype giniimpurity(subset,classes)

@show extrema(filter(isfinite,substrate_index))

substrate_index_int = map(x -> (isfinite(x) ? round(Int,x) : 0), substrate_index)

@show extrema(filter(isfinite,substrate_index_int))

classes = 0:nsubstrate

i = 50:70
j = 50:70
ibox = 10
ibox = 50
ibox = 20

function substrate_gini_impurity(substrate_index_int,classes,ibox)
    sz = size(substrate_index_int)

    gimpurity = zeros(sz)

    for j = 1:size(substrate_index_int,2)
        @show j
        for i = 1:size(substrate_index_int,1)
            si = max(i-ibox,1):min(i+ibox,sz[1])
            sj = max(j-ibox,1):min(j+ibox,sz[2])
            gimpurity[i,j] = giniimpurity(view(substrate_index_int,si,sj),classes)
        end
    end
    return gimpurity
end

gimpurity = @time substrate_gini_impurity(substrate_index_int,classes,ibox)
size(substrate_index)

#gimpurity_ref = copy(gimpurity)
clf();pcolormesh(s_lon,s_lat,gimpurity')
colorbar()
OceanPlot.plot_coastline()

typeof(gimpurity)

gimpurity2 = replace(gimpurity,NaN => missing)
Dataset(joinpath(datadir,"substrate_gini_impurity.nc"),"c") do ds
    defVar(ds,"lon",s_lon,("lon",))
    defVar(ds,"lat",s_lat,("lat",))
    defVar(ds,"substrate_gini_impurity",gimpurity,("lon","lat"))
end



maxlen = 100e3
minlen = 20e3

len = minlen .+ (maxlen - minlen) * (1 .- gimpurity)
clf();pcolormesh(s_lon,s_lat,len')
colorbar()
OceanPlot.plot_coastline()
title("Sample correlation length (m)")
savefig("../figures/corr_len_test.png")
