"""Download ERA5 pressure-level monthly means for 2024 (all 12 months).

Fetches only the variables tcpyPI consumes, subset to the North Atlantic region,
so the download stays small (~100 MB rather than several GB). Used to build the
2024 sample (data/sample_data.nc) and the single-column demo subset.

Requires a CDS API key in ~/.cdsapirc (https://cds.climate.copernicus.eu/how-to-api).
Writes: data/era5_pl_monthlymeans_2024.nc  (run from the repo's data/ directory).
"""

import cdsapi

dataset = "reanalysis-era5-pressure-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "temperature",  # t  -> TC (after K->C)
        "specific_humidity",  # q  -> R  (g/kg)
        "u_component_of_wind",  # u  (shear, hodograph)
        "v_component_of_wind",  # v  (shear, hodograph)
        "relative_humidity",  # r  (GPI en04 midlevel RH)
        "vorticity",  # vo (GPI absolute vorticity)
    ],
    "pressure_level": [
        "50",
        "70",
        "100",
        "125",
        "150",
        "175",
        "200",
        "225",
        "250",
        "300",
        "350",
        "400",
        "450",
        "500",
        "550",
        "600",
        "650",
        "700",
        "750",
        "775",
        "800",
        "825",
        "850",
        "875",
        "900",
        "925",
        "950",
        "975",
        "1000",
    ],
    "year": ["2024"],
    "month": [
        "01",
        "02",
        "03",
        "04",
        "05",
        "06",
        "07",
        "08",
        "09",
        "10",
        "11",
        "12",
    ],
    "time": ["00:00"],
    # North Atlantic region [North, West, South, East] (lon in -180..180):
    # covers the main development region + Gulf of Mexico (Milton).
    "area": [50, -110, 0, -40],
    "data_format": "netcdf",
    "download_format": "unarchived",
}

client = cdsapi.Client()
client.retrieve(dataset, request, "era5_pl_monthlymeans_2024.nc")
