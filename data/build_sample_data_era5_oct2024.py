"""Build the pyPI 2024 sample datasets from ERA5 monthly means.

Derives a small North Atlantic regional subset suitable for running `run_sample.py`
and the seasonal analyses, plus a single-column demo subset for the ventilation-index
and GPI notebooks. Handles any number of monthly-mean time steps, so the same script
builds a 12-month seasonal sample or a single-month sample.

Inputs (expected in `data/`, produced by the download scripts):
- `era5_pl_monthlymeans_2024.nc` (pressure-level monthly means; from download_era5.py)
- `era5_sl_monthlymeans_2024.nc` (single-level monthly means; from download_era5_sl.py)
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

from tcpyPI.utilities import convert_lon_to360


def _subset_latlon(ds: xr.Dataset, lat_min: float, lat_max: float, lon_min: float, lon_max: float) -> xr.Dataset:
    lat = ds["lat"]
    lon = ds["lon"]

    ds = ds.where((lat >= lat_min) & (lat <= lat_max), drop=True)

    if lon_min <= lon_max:
        ds = ds.where((lon >= lon_min) & (lon <= lon_max), drop=True)
        return ds

    ds = ds.where((lon >= lon_min) | (lon <= lon_max), drop=True)
    return ds


def build_sample_data(
    pl_path: str | Path = Path("data") / "era5_pl_monthlymeans_2024.nc",
    sl_path: str | Path = Path("data") / "era5_sl_monthlymeans_2024.nc",
    out_path: str | Path = Path("data") / "sample_data.nc",
    *,
    lat_min: float = 0.0,
    lat_max: float = 50.0,
    lon_min: float = 250.0,
    lon_max: float = 320.0,
    subsample_factor: int = 8,
    land_frac_max: float = 0.5,
    p_min: float = 50.0,
    p_max: float = 1000.0,
) -> Path:
    """Create `sample_data.nc` in the pyPI format (C/hPa/g/kg; dims `month,p,lat,lon`).

    Handles any number of monthly-mean time steps: the ``valid_time`` axis is mapped
    to an integer ``month`` coordinate (1-12), so a 12-month 2024 download yields a
    full seasonal cycle and a single-month download yields a one-element ``month``
    axis. Includes the wind/humidity/vorticity fields (``u, v, r, vo``) needed for the
    seasonal GPI/ventilation-index analyses alongside the core PI inputs
    (``sst, msl, t, q``).

    Spatial resolution is reduced by strided subsampling (decimation, no
    interpolation) with ``subsample_factor`` (8 -> ~2 deg from ERA5's 0.25 deg); the
    full native field is always regenerable with ``subsample_factor=1``. SST is masked
    to ocean where the ERA5 land fraction is >= ``land_frac_max`` (default 0.5), so PI
    is not computed over (partly) land columns.
    """

    pl_path = Path(pl_path)
    sl_path = Path(sl_path)
    out_path = Path(out_path)

    pl = xr.open_dataset(pl_path)[["t", "q", "u", "v", "r", "vo"]]
    pl = pl.rename({"pressure_level": "p", "latitude": "lat", "longitude": "lon"})
    pl = pl.drop_vars(["number", "expver"], errors="ignore")
    # Normalize longitude to 0-360 (a regional `area` download returns -180..180),
    # so the lon subset below and the output convention are convention-agnostic.
    pl = pl.assign_coords(lon=convert_lon_to360(pl["lon"])).sortby("lon")
    # Sort to descending pressure first so the level subset is order-agnostic
    # (slice(p_max, p_min) returns empty on an ascending-pressure source otherwise).
    pl = pl.sortby("p", ascending=False).sel(p=slice(p_max, p_min))
    pl = _subset_latlon(pl, lat_min, lat_max, lon_min, lon_max)

    sl = xr.open_dataset(sl_path)[["sst", "msl", "lsm"]]
    sl = sl.rename({"latitude": "lat", "longitude": "lon"})
    sl = sl.drop_vars(["number", "expver"], errors="ignore")
    sl = sl.assign_coords(lon=convert_lon_to360(sl["lon"])).sortby("lon")
    sl = _subset_latlon(sl, lat_min, lat_max, lon_min, lon_max)

    # Reduce resolution by strided subsampling (decimation) — keeps true native ERA5
    # gridpoint values (no interpolation/averaging). Use xESMF if conservative
    # regridding is ever needed instead.
    if subsample_factor > 1:
        sub = {"lat": slice(None, None, subsample_factor), "lon": slice(None, None, subsample_factor)}
        pl = pl.isel(sub)
        sl = sl.isel(sub)

    # Map the monthly-means time axis to an integer month coordinate (1-12).
    months_pl = pl["valid_time"].dt.month.values.astype("int64")
    months_sl = sl["valid_time"].dt.month.values.astype("int64")
    pl = pl.rename({"valid_time": "month"}).assign_coords(month=("month", months_pl))
    sl = sl.rename({"valid_time": "month"}).assign_coords(month=("month", months_sl))

    # Mask SST to ocean: drop cells that are >= land_frac_max land, so PI is not
    # computed over (partly) land columns (which yield contaminated, inflated PI).
    sst_c = (sl["sst"] - 273.15).where(sl["lsm"] < land_frac_max)

    out = xr.Dataset(
        data_vars={
            # lsm is time-invariant; store it 2-D (lat, lon).
            "lsm": sl["lsm"].isel(month=0, drop=True),
            "sst": sst_c,
            "msl": sl["msl"] / 100.0,
            "t": pl["t"] - 273.15,
            "q": pl["q"] * 1000.0,
            "u": pl["u"],
            "v": pl["v"],
            "r": pl["r"],
            "vo": pl["vo"],
        },
        attrs={
            "title": "pyPI sample dataset (ERA5 monthly means, 2024)",
            "source": f"{pl_path.name} + {sl_path.name}",
            "history": f"{datetime.now(timezone.utc).isoformat()}: created by data/build_sample_data_era5_oct2024.py",
            "notes": (
                "Core PI inputs in pyPI units: SST/T in C, MSL in hPa, q in g/kg "
                "(specific humidity used directly as the mixing-ratio input, following "
                "the original pyPI sample convention). u/v [m/s], r [%], vo [s^-1] are "
                "included for the seasonal GPI/ventilation-index analyses. Native ERA5 "
                "0.25deg subsampled (decimated) to ~2deg; SST masked to ocean "
                "(land fraction < 0.5)."
            ),
        },
    )

    out["month"].attrs.update({"standard_name": "Month", "units": "month number (1-12)"})
    out["p"].attrs.update({"standard_name": "Atmospheric Pressure", "units": "hPa"})
    out["lat"].attrs.update({"standard_name": "Latitude", "units": "degrees_north"})
    out["lon"].attrs.update({"standard_name": "Longitude", "units": "degrees_east"})

    out["lsm"].attrs.update({
        "standard_name": "land_area_fraction",
        "long_name": "Land fraction (native ERA5 land-sea mask at the subsampled points)",
        "units": "1",
    })
    out["sst"].attrs.update({"standard_name": "Sea Surface Temperature", "units": "degrees C"})
    out["msl"].attrs.update({"standard_name": "Mean Sea Level Pressure", "units": "hPa"})
    out["t"].attrs.update({"standard_name": "Atmospheric Temperature", "units": "degrees C"})
    out["q"].attrs.update({"standard_name": "Specific Humidity", "units": "g/kg"})
    out["u"].attrs.update({"standard_name": "Eastward Wind", "units": "m/s"})
    out["v"].attrs.update({"standard_name": "Northward Wind", "units": "m/s"})
    out["r"].attrs.update({"standard_name": "Relative Humidity", "units": "%"})
    out["vo"].attrs.update({"standard_name": "Relative Vorticity", "units": "s-1"})

    encoding = {var: {"zlib": True, "complevel": 4} for var in out.data_vars}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(out_path, encoding=encoding, engine="h5netcdf")
    return out_path


def build_era5_demo_subset(
    pl_path: str | Path = Path("data") / "era5_pl_monthlymeans_2024.nc",
    sl_path: str | Path = Path("data") / "era5_sl_monthlymeans_2024.nc",
    out_path: str | Path = Path("data") / "era5_demo_subset.nc",
    *,
    milton_lat: float = 21.8,
    milton_lon: float = 269.1,
    milton_time: str = "2024-10-07 20:00 UTC",
    p_min: float = 50.0,
    p_max: float = 1000.0,
) -> Path:
    """Build a small, tracked single-column ERA5 demo environment.

    Extracts the nearest ERA5 October-2024 monthly-mean column at the Hurricane
    Milton (2024) peak-intensity reference point (default 21.8 N, 269.1 E, from
    IBTrACS) for the ventilation-index and GPI demo notebooks. Kept small and
    committed to the repo so those notebooks run on a fresh clone without the full
    ERA5 download or IBTrACS.

    Variables are left in native ERA5 units (``t`` [K], ``q`` [kg/kg], ``u``/``v``
    [m/s], ``r`` [%], ``vo`` [s^-1], ``sst`` [K], ``msl`` [Pa]); the notebooks
    convert as needed. The Milton reference lat/lon/time are stored in the dataset
    attributes so notebooks need not re-read IBTrACS.
    """
    pl_path = Path(pl_path)
    sl_path = Path(sl_path)
    out_path = Path(out_path)

    # Select the October time step (Milton era); works for a 12-month file or a
    # single-month October file.
    pl_full = xr.open_dataset(pl_path)[["t", "q", "u", "v", "r", "vo"]]
    # Normalize longitude to 0-360 (a regional `area` download returns -180..180).
    pl_full = pl_full.assign_coords(longitude=convert_lon_to360(pl_full["longitude"])).sortby("longitude")
    oct_idx = int(np.argmax(pl_full["valid_time"].dt.month.values == 10))
    pl = pl_full.isel(valid_time=oct_idx, drop=True)
    pl = pl.drop_vars(["number", "expver"], errors="ignore")
    pl = pl.sortby("pressure_level", ascending=False).sel(
        pressure_level=slice(p_max, p_min)
    )

    sl_full = xr.open_dataset(sl_path)[["sst", "msl"]]
    sl_full = sl_full.assign_coords(longitude=convert_lon_to360(sl_full["longitude"])).sortby("longitude")
    sl_oct = int(np.argmax(sl_full["valid_time"].dt.month.values == 10))
    sl = sl_full.isel(valid_time=sl_oct, drop=True)
    sl = sl.drop_vars(["number", "expver"], errors="ignore")

    pl_pt = pl.sel(latitude=milton_lat, longitude=milton_lon, method="nearest")
    sl_pt = sl.sel(latitude=milton_lat, longitude=milton_lon, method="nearest")

    out = xr.merge([pl_pt, sl_pt])
    out.attrs.update(
        {
            "title": "pyPI ERA5 demo subset (single column; Hurricane Milton 2024 reference)",
            "source": f"{pl_path.name} + {sl_path.name}",
            "reference": "Location near the Hurricane Milton (2024) IBTrACS peak-intensity point",
            "milton_lat": float(pl_pt["latitude"].values),
            "milton_lon": float(pl_pt["longitude"].values),
            "milton_lon_convention": "degrees_east (0-360)",
            "milton_time": milton_time,
            "units_note": (
                "Native ERA5 units: t[K], q[kg/kg], u/v[m/s], r[%], vo[s^-1], "
                "sst[K], msl[Pa]."
            ),
            "caveat": (
                "This is an October-2024 MONTHLY-MEAN (climatological) column - a valid "
                "basis for climatological PI/GPI/VI diagnostics across many storms or "
                "seasons. It is NOT the pre-storm environment of this individual storm: "
                "to characterize the environment a specific storm developed in, use "
                "sub-monthly fields 3-5 days AHEAD of passage so the storm does not "
                "contaminate its own environment. Milton's location is used here only as "
                "a concrete example."
            ),
            "notes": (
                "Single-column environment for demonstrating the ventilation-index and "
                "GPI calculations in the demo notebooks; see "
                "data/build_sample_data_era5_oct2024.py."
            ),
        }
    )

    encoding = {var: {"zlib": True, "complevel": 4} for var in out.data_vars}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(out_path, encoding=encoding, engine="h5netcdf")
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the pyPI 2024 sample and demo subset from ERA5 monthly means."
    )
    parser.add_argument("--pl", default="data/era5_pl_monthlymeans_2024.nc", help="Pressure-level ERA5 NetCDF")
    parser.add_argument("--sl", default="data/era5_sl_monthlymeans_2024.nc", help="Single-level ERA5 NetCDF")
    parser.add_argument("--out", default="data/sample_data.nc", help="Output NetCDF path")
    parser.add_argument(
        "--demo-out",
        default="data/era5_demo_subset.nc",
        help="Output path for the small single-column VI/GPI demo subset (set empty to skip)",
    )
    parser.add_argument("--lat-min", type=float, default=0.0)
    parser.add_argument("--lat-max", type=float, default=50.0)
    parser.add_argument("--lon-min", type=float, default=250.0)
    parser.add_argument("--lon-max", type=float, default=320.0)
    parser.add_argument("--subsample", type=int, default=8, help="Strided decimation factor for lat/lon (8 -> ~2 deg from 0.25; 1 = full native)")
    parser.add_argument("--land-frac-max", type=float, default=0.5, help="Mask SST where land fraction >= this (ocean-only PI)")
    parser.add_argument("--p-min", type=float, default=50.0)
    parser.add_argument("--p-max", type=float, default=1000.0)
    args = parser.parse_args()

    out_path = build_sample_data(
        args.pl,
        args.sl,
        args.out,
        lat_min=args.lat_min,
        lat_max=args.lat_max,
        lon_min=args.lon_min,
        lon_max=args.lon_max,
        subsample_factor=args.subsample,
        land_frac_max=args.land_frac_max,
        p_min=args.p_min,
        p_max=args.p_max,
    )
    print(f"Wrote {out_path}")

    if args.demo_out:
        demo_path = build_era5_demo_subset(
            args.pl,
            args.sl,
            args.demo_out,
            p_min=args.p_min,
            p_max=args.p_max,
        )
        print(f"Wrote {demo_path}")


if __name__ == "__main__":
    main()

