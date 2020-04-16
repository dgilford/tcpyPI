# This pyPI script computes PI and associated analyses over the entire sample dataset
# which is from 2004, MERRA2.
#
# Created by Daniel Gilford, PhD (daniel.gilford@rutgers.edu)
# Many thanks to Daniel Rothenberg for his assitance optimizing pyPI
#
# Last updated 4/16/2020
#

# setup
import xarray as xr
import pickle

# load in pyPI modules
from pi import pi

# define the sample data locations
datdir='./data/'
_FN=datdir+'sample_data.nc'
_mdrF=datdir+'mdr.pk1' 
    

def run_sample_dataset(fn, dim='p',CKCD=0.9):
    """ This function calculates PI over the sample dataset using xarray """
    
    # open the sample data file
    ds = xr.open_dataset(fn)
    # calculate PI over the whole data set using the xarray universal function
    result = xr.apply_ufunc(
        pi,
        ds['sst'], ds['msl'], ds['p'], ds['t'], ds['q'],
        kwargs=dict(CKCD=CKCD, ascent_flag=0, diss_flag=1, miss_handle=1),
        input_core_dims=[
            [], [], ['p', ], ['p', ], ['p', ],
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
        'q': ds.q,
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

def run_sample_analyses(ds,_mdrF,CKCD=0.9):
    """ This function performs PI analyses over the sample dataset using xarray """

    # load the basins dictionary
    basins = pickle.load( open( _mdrF, "rb" ) )
    
    # import functions for analyses
    from utilities import pi_effiency
    from utilities import pi_diseq_resid
    from utilities import decompose_pi
    
    # calculate PI analyses over the whole data set using the xarray universal function
    efficiency = xr.apply_ufunc(
        pi_effiency,
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
        decompose_pi,
        ds['vmax'], ds['sst']+273.15, ds['t0'],
        kwargs=dict(CKCD=CKCD),
        input_core_dims=[
            [], [], [],
        ],
        output_core_dims=[
            [], [], [], [],
        ],
        vectorize=True
    )

    lnpi, lneff, lndiseq, lnCKCD = result
    
    out_ds = xr.Dataset({
                'eff': efficiency, 
                'diseq': diseq,
                'lnpi': lnpi,
                'lneff': lneff,
                'lndiseq': lndiseq,
                'lnCKCD': lnCKCD[0,0,0]
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
    ds.to_netcdf(datdir+'raw_sample_output.nc')
    print('...PI computation complete and saved\n')
    
    # Perform PI analyses over the whole dataset
    print('Performing PI analyses...')
    ds2 = run_sample_analyses(ds,_mdrF,CKCD=0.9)
    
    # merge the arrays and save the output
    ds3=ds.merge(ds2)
    ds3.to_netcdf(datdir+'full_sample_output.nc')
    del ds, ds2
    print('...PI analyses complete and saved')
    print('pyPI sample run finished!')