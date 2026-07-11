"""Download ERA5 single-level monthly means for 2024 (all 12 months).

Fetches only the fields tcpyPI consumes (SST, MSL, land-sea mask), subset to the
North Atlantic region. Requesting only these analysis-stream variables keeps CDS
from returning a multi-file zip, so the result is a single netCDF.

Requires a CDS API key in ~/.cdsapirc (https://cds.climate.copernicus.eu/how-to-api).
Writes: data/era5_sl_monthlymeans_2024.nc  (run from the repo's data/ directory).
"""

import cdsapi

dataset = "reanalysis-era5-single-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "sea_surface_temperature",   # sst -> SSTC (after K->C)
        "mean_sea_level_pressure",   # msl -> MSL  (after Pa->hPa)
        "land_sea_mask",             # lsm (land fraction after coarsening)
    ],
    "year": ["2024"],
    "month": [
        "01", "02", "03", "04", "05", "06",
        "07", "08", "09", "10", "11", "12",
    ],
    "time": ["00:00"],
    # North Atlantic region [North, West, South, East] (lon in -180..180).
    "area": [50, -110, 0, -40],
    "data_format": "netcdf",
    "download_format": "unarchived",
}

client = cdsapi.Client()
client.retrieve(dataset, request, "era5_sl_monthlymeans_2024.nc")
