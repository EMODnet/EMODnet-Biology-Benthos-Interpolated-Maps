# EMODnet-Biology-Benthos-Interpolated-Maps

This directory provides the code to
* read the benthos data files provided by P. Hermann and
* perform spatial interpolations using `DIVAnd` software tool.

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
* **product** - Output product files: `netCDF` containing the gridded, probability fields and the corresponding figures in `PNG` format.
* **scripts** - Reusable code: functions employed in the Jupyter notebooks.

## Data

All the data records are stored in the file `specs4Diva.csv`.       

The data is structured like this:
```bash
eventNummer,eventDate,decimalLongitude,decimalLatitude,Abra_alba,Amphiura_filiformis,Diplocirrus_glaucus,Amphiura_chiajei,Abludomelita_obtusata,Bicellariella_ciliata,Megaluropus_agilis,Sigalion_mathildae,Callianassa_subterranea,Acrocnida_brachiata,Littorina_littorea,Aequipecten_opercularis,Bathyporeia_tenuipes,Lumbrineriopsis_paradoxa,Spio_armata,Schizomavella_(Schizomavella)_auriculata,Amphiura_(Ophiopeltis)_securigera,Diplocirrus_stopbowitzi
1,2000-01-01T00:00:00Z,8.064666670000005,54.76183333000002,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
2,2001-01-01T00:00:00Z,8.420666670000005,54.28333333000003,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
3,2000-01-01T00:00:00Z,8.01166667,54.79633333000006,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
4,2000-01-01T00:00:00Z,8.09766667,54.731333330000034,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
...
```
so we have:
* the event number,
* the time of the observation,
* the coordinates and
* a `TRUE` or `FALSE` value for different species.





## Analysis


## Citation

Please cite this product as:
*A. Barth, P. Hermann and C. Troupin (2020). Probability maps
for different benthos species in the North Sea.*
