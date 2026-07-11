# This pyPI script computes PI and associated analyses over the sample dataset.
#
# Created by Daniel Gilford, PhD (daniel.gilford@rutgers.edu)
# Many thanks to Daniel Rothenberg for his assitance optimizing pyPI
#
# Last updated 8/14/2020
#

# setup
import sys
from pathlib import Path

import xarray as xr

# Prefer the local `src/` checkout when running from a repo clone.
_REPO_ROOT = Path(__file__).resolve().parent
_SRC_DIR = _REPO_ROOT / "src"
if _SRC_DIR.is_dir():
    sys.path.insert(0, str(_SRC_DIR))

# load in pyPI modules
from tcpyPI import pi, pi_log_decomposition
from tcpyPI.utilities import pi_efficiency, pi_diseq_resid

# define the sample data locations
datdir = "./data/"
_FN = datdir + "sample_data.nc"
    

def run_sample_dataset(fn, dim='p',CKCD=0.9,outflow_source="cape_star"):
    """Calculate potential intensity over the sample dataset with xarray.

    Parameters
    ----------
    fn : str
        Path to the input dataset (e.g., `data/sample_data.nc`).
    dim : str, default="p"
        Name of the vertical pressure coordinate in `fn`.
    CKCD : float, default=0.9
        Ratio of exchange coefficients (Ck/Cd).
    outflow_source : {"cape_star", "cape_env"}, default="cape_star"
        Which CAPE calculation supplies the outflow temperature and level.

    Returns
    -------
    xarray.Dataset
        Dataset containing PI outputs (`vmax`, `pmin`, `ifl`, `t0`, `otl`) and the
        state variables used in the calculation.
    """
    
    # open the sample data file
    ds = xr.open_dataset(fn)
    # calculate PI over the whole data set using the xarray universal function
    result = xr.apply_ufunc(
        pi,
        ds['sst'], ds['msl'], ds[dim], ds['t'], ds['r'],
        kwargs=dict(
            CKCD=CKCD,
            ascent_flag=0,
            diss_flag=1,
            ptop=50,
            miss_handle=1,
            outflow_source=outflow_source,
        ),
        input_core_dims=[
            [], [], [dim, ], [dim, ], [dim, ],
        ],
        output_core_dims=[
            [], [], [], [], []
        ],
        vectorize=True
    )

    # store the result in an xarray data structure
    vmax, pmin, ifl, t0, otl = result
    out_ds=xr.Dataset({
        'vmax': vmax, 
        'pmin': pmin,
        'ifl': ifl,
        't0': t0,
        'otl': otl,
        # merge the state data into the same data structure
        'sst': ds.sst,
        't': ds.t,
        'r': ds.r,
        'msl': ds.msl,
        'lsm': ds.lsm,
        })
    
    # add names and units to the structure
    out_ds.vmax.attrs['standard_name'],out_ds.vmax.attrs['units']='Maximum Potential Intensity','m/s'
    out_ds.pmin.attrs['standard_name'],out_ds.pmin.attrs['units']='Minimum Central Pressure','hPa'
    out_ds.ifl.attrs['standard_name']='pyPI Flag'
    out_ds.t0.attrs['standard_name'],out_ds.t0.attrs['units']='Outflow Temperature','K'
    out_ds.otl.attrs['standard_name'],out_ds.otl.attrs['units']='Outflow Temperature Level','hPa'

    # return the output from pi.py as an xarray data structure
    return out_ds

def run_sample_analyses(ds,CKCD=0.9):
    """Compute analysis diagnostics (including log decomposition) on PI outputs.

    Parameters
    ----------
    ds : xarray.Dataset
        Output from :func:`run_sample_dataset`.
    CKCD : float, default=0.9
        Ratio of exchange coefficients (Ck/Cd).

    Returns
    -------
    xarray.Dataset
        Dataset of diagnostics: `eff`, `diseq`, `lnpi`, `lneff`, `lndiseq`, `lnCKCD`.
    """
    
    # calculate PI analyses over the whole data set using the xarray universal function
    efficiency = xr.apply_ufunc(
        pi_efficiency,
        ds['sst']+273.15, ds['t0'],
        input_core_dims=[
            [], [],
        ],
        output_core_dims=[
            [],
        ],
        vectorize=True
    )
    
    diseq = xr.apply_ufunc(
        pi_diseq_resid,
        ds['vmax'], ds['sst']+273.15, ds['t0'],
        kwargs=dict(CKCD=CKCD),
        input_core_dims=[
            [], [], [],
        ],
        output_core_dims=[
            [],
        ],
        vectorize=True
    )
    
    result = xr.apply_ufunc(
        pi_log_decomposition,
        ds['vmax'], ds['sst']+273.15, ds['t0'],
        kwargs=dict(CKCD=CKCD, sst_units="K"),
        input_core_dims=[
            [], [], [],
        ],
        output_core_dims=[
            [], [], [], [],
        ],
        vectorize=True
    )

    lnpi, lneff, lndiseq, lnCKCD = result
    
    # `lnCKCD` is constant in space/time, but `pi_log_decomposition` broadcasts it to
    # the input shape. Collapse it to a scalar so the output stays compact and
    # robust to input dimensionality.
    lnCKCD_scalar = lnCKCD
    for dim in list(getattr(lnCKCD_scalar, "dims", ())):
        lnCKCD_scalar = lnCKCD_scalar.isel({dim: 0})

    out_ds = xr.Dataset({
                'eff': efficiency, 
                'diseq': diseq,
                'lnpi': lnpi,
                'lneff': lneff,
                'lndiseq': lndiseq,
                'lnCKCD': lnCKCD_scalar
            })
    
    # add names and units (where applicable)
    out_ds.eff.attrs['standard_name'],out_ds.eff.attrs['units']='Tropical Cyclone Efficiency','unitless fraction'
    out_ds.diseq.attrs['standard_name'],out_ds.diseq.attrs['units']='Thermodynamic Disequilibrium','J/kg'
    out_ds.lnpi.attrs['standard_name']='Natural log(Potential Intensity)'
    out_ds.lneff.attrs['standard_name']='Natural log(Tropical Cyclone Efficiency)'
    out_ds.lndiseq.attrs['standard_name']='Natural log(Thermodynamic Disequilibrium)'
    out_ds.lnCKCD.attrs['standard_name'],out_ds.lnCKCD.attrs['units']='Natural log(Ck/CD)','unitless constant'

    # return the output from pi.py as an xarray data structure
    return out_ds
    
    

if __name__ == "__main__":

    # Execute PI analysis over the whole dataset and save the output
    print('Beginning PI computations...')
    ds = run_sample_dataset(_FN)
    ds.to_netcdf(datdir+'raw_sample_output.nc', engine="h5netcdf")
    print('...PI computation complete and saved\n')
    
    # Perform PI analyses over the whole dataset
    print('Performing PI analyses...')
    ds2 = run_sample_analyses(ds,CKCD=0.9)
    
    # merge the arrays and save the output
    ds3=ds.merge(ds2)
    ds3.to_netcdf(datdir+'full_sample_output.nc', engine="h5netcdf")
    del ds, ds2
    print('...PI analyses complete and saved')
    print('pyPI sample run finished!')
