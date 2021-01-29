# EMODnet-Biology-Benthos-Interpolated-Maps

`DIVA` (Data Interpolating Variational Analysis) and `DIVAnd` (DIVA in __n__ dimensions) are software tools primarily designed to generate gridded maps of continuous variables such as sea water temperature, salinity or oxygen concentration. The main advantages over other interpolation or analysis method are:
* coastlines and physical boundaries are taken into account by the method.
* large datasets (million of data points) can be ingested and processed by the tool.

DIVAnd is a multi-dimensional generalization ([Barth et al., 2014](https://dx.doi.org/10.5194/gmd-7-225-2014)), written in the [Julia language](https://julialang.org/), with a new mathematical formulation with respect to the previous DIVA code.

This directory provides the codes and tools to
- process the presence/absence data relative to different _benthos_ species, and
- generate gridded field using `DIVAnd`.

## DIVAnd for presence/absence data

The dataset to be processed contained only a binary information: presence or absence, and the objective is to derive a map showing the probability to encounter a given species, as illustrated below for _Aequipecten opercularis_.

![presense/absence](product/figures/1-UniformL/data/Aequipecten_opercularis_data.jpg)


## Directory structure

```
EMODnet-Biology-Benthos-Interpolated-Maps/
├── analysis
├── data/
├── docs/
├── product/
│   ├── figures/
│   └── netCDF/
└── scripts/
```

* **analysis** - Jupyter notebooks used to perform the data analysis, create the figures and the `netCDF` files.
* **data** - input data files.
* **docs** - Rendered reports
* **product** - Output product files: `netCDF` containing the gridded, probability fields and the corresponding figures in `PNG` or `JPG` format.
* **scripts** - Reusable code: functions employed in the Jupyter notebooks.

## Data

The data records stored in the file `specs4Diva.csv` consists of 18 species of benthos that were used for testing purposes.     
The main data file is `spe.csv`, it contains all sampling events as rows, and the most frequent species (occurrence >200) as columns. The first columns specify sampling event and location. All columns are coded as ‘pa’ (for presence/absence) followed by the AphiaID of the species. In the separate file specieslist.csv you find the names of the species corresponding to their AphiaID. You also find the total number of events where the species has been found.

A species can have value 1 (present), 0 (looked for but absent) and NA (not looked for). The NA values must be disregarded when making maps. They do not contain information about either presence or absence. There are NAs in the file because some datasets only look for a subset of species. 

The data file is structured like this:
```bash
eventNummer,eventDate,decimalLongitude,decimalLatitude,Abra_alba,Amphiura_filiformis,Diplocirrus_glaucus,Amphiura_chiajei,Abludomelita_obtusata,Bicellariella_ciliata,Megaluropus_agilis,Sigalion_mathildae,Callianassa_subterranea,Acrocnida_brachiata,Littorina_littorea,Aequipecten_opercularis,Bathyporeia_tenuipes,Lumbrineriopsis_paradoxa,Spio_armata,Schizomavella_(Schizomavella)_auriculata,Amphiura_(Ophiopeltis)_securigera,Diplocirrus_stopbowitzi
1,2000-01-01T00:00:00Z,8.064666670000005,54.76183333000002,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
2,2001-01-01T00:00:00Z,8.420666670000005,54.28333333000003,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
3,2000-01-01T00:00:00Z,8.01166667,54.79633333000006,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
4,2000-01-01T00:00:00Z,8.09766667,54.731333330000034,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
...
```
so we have:
* the event number (not used in the present analysis),
* the time of the observation (not used),
* the coordinates (WGS84 projection EPGS 4326) and
* a `TRUE` (presence) or `FALSE` (absence)  value for the 18 different species.

The information is read by the function
`read_coords_species` defined in the script  [`BenthosInterp.jl`](scripts/BenthosInterp.jl).

## Analysis

This directory contains the notebooks for the preparation and analysis of the data, as well as the generation of the figures.

* `interp_presence_absence_benthos.ipynb`:
1. prepare the domain and the bathymetry,
2. read the coordinates of presence and absence for each species,
3. compute the probability map and the associated error field and
4. write the results in a netCDF file (one per species), in the directory `product/netCDF`

* `plot_results_map.ipynb`: notebook in Python to create the figures using the _ETRS89 Lambert Azimuthal Equal Area_ coordinate reference system of 2001 (`EPGS 3035`).     
The figures are stored in `product/figures` and consist of
1. The interpolated probability map
2. The relative error field
3. The interpolated probability masked by the error field

![variableL](product/figures/1-UniformL/)
![variableL](product/figures/1-UniformL/)


### Products

Two types of analysis are performed:
1. Using a uniform correlation length all over the domain (0.1°) (directory `1-UniformL`)
2. Using a spatially variable correlation length, derives from the substrates (directory
	`2-VariableL`).

![variableL](product/figures/variableL.jpg)

### Code

The code, written in Julia, is distributed through GitHub:
https://github.com/EMODnet/EMODnet-Biology-Benthos-Interpolated-Maps
The `analysis` directory stores different jupyter-notebooks describing the different analysis steps.

## Citation

Please cite this product as:
*A. Barth, P. Hermann and C. Troupin (2020). Probability maps
for different benthos species in the North Sea. Integrated data products created under the European Marine Observation Data Network (EMODnet) Biology project (EASME/EMFF/2017/1.3.1.2/02/SI2.789013), funded by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund.*

The code is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
