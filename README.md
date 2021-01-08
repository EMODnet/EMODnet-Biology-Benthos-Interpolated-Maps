# EMODnet-Biology-Benthos-Interpolated-Maps

`DIVA` (Data Interpolating Variational Analysis) and `DIVAnd`ß (DIVA in n dimensions) are software tools primarily designed to generate gridded maps of continuous variables such as sea water temperature, salinity or oxygen concentration. The main advantages over other interpolation or analysis method are:
* coastlines and physical boundaries are taken into account by the method.
* large datasets (million of data points) can be ingested and processed by the tool.

DIVAnd is a multi-dimensional generalization ([Barth et al., 2014](https://dx.doi.org/10.5194/gmd-7-225-2014)), written in the [Julia language](https://julialang.org/), with a new mathematical formulation with respect to the previous DIVA code.

This directory provides the codes and tools to
- process the presence/absence data relative to different benthos species, and
- generat gridded field using `DIVAnd`.

## DIVAnd for presence/absence data

The dataset to be processed contained only a binary information: presence or absence, and the objective is to derive a map showing the probability to encounter a given species, as illustrated below for _Aequipecten opercularis_.

![presense/absence](product/figures/1-UniformL/data/Aequipecten_opercularis_data.jpg)

* [Heatmap](analysis/interp_presence_absence_benthos.ipynb):


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

All the data records are stored in the file `specs4Diva.csv` and contain 18 species of benthos.       

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


## Citation

Please cite this product as:
*A. Barth, P. Hermann and C. Troupin (2020). Probability maps
for different benthos species in the North Sea.*
