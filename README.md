# tcpyPI: Potential Intensity Calculations in Python
<p align="center">
<img src="./figures/readme_image.png" alt="" width="720" height="480">
</p>

tcpyPI, 'pyPI' for short, is a set of scripts and notebooks that compute and validate tropical cyclone (TC) potential intensity (PI) calculations in Python.
It is a fully documented and improved port of the [Bister and Emanuel 2002](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2001JD000776) algorithm (hereafter BE02) which was originally written in FORTRAN---and then MATLAB---by Prof. Kerry Emanuel (MIT).
Kerry's original MATLAB code (**pcmin.m**) is found at:

* [http://texmex.mit.edu/pub/emanuel/TCMAX](http://texmex.mit.edu/pub/emanuel/TCMAX)

The goals in developing and maintaining pyPI are to:
* supply a freely available validated Python potential intensity calculator,
* carefully document the BE02 algorithm and its Python implementation, and to
* demonstrate and encourage the use of potential intensity theory in tropical cyclone climatology analysis.

If you have any questions, comments, or feedback, please [contact the developer](mailto:dgilford@climatecentral.org) or open an [Issue](https://github.com/dgilford/tcpyPI/issues) in the repository. A paper detailing pyPI is published [at Geoscientific Model Development](https://gmd.copernicus.org/articles/14/2351/2021/gmd-14-2351-2021.pdf).

## Citation
pyPI was developed by [Daniel Gilford](https://github.com/dgilford) and is archived on Zenodo. A machine-readable [`CITATION.cff`](CITATION.cff) is included.

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.3756005-1682D4.svg)](https://doi.org/10.5281/zenodo.3756005)

If you use pyPI in your work, please cite **both** the paper and the software. The Zenodo **concept DOI** below always resolves to the latest version (cite a specific version DOI if you need it, e.g. `10.5281/zenodo.21301851` for v1.4.2):

> Gilford, D. M.: pyPI (v1.3): Tropical Cyclone Potential Intensity Calculations in Python, Geosci. Model Dev., 14, 2351–2369, https://doi.org/10.5194/gmd-14-2351-2021, 2021.

and

> Gilford, D. M.: tcpyPI (pyPI): Tropical Cyclone Potential Intensity Calculations in Python. Zenodo. https://doi.org/10.5281/zenodo.3756005


## Full pyPI Description

Please read [pyPI_Users_Guide_v1.3.pdf](pyPI_Users_Guide_v1.3.pdf) for a full overview and details on pyPI.
The description includes the pyPI background, a PI computation derivation, validation against the commonly-used MATLAB algorithm (pcmin), and a set of sample analyses.

> **Note:** the User's Guide documents v1.3 and predates the newer options and diagnostics (`outflow_source`, the log-decomposition API, and the ventilation, genesis-potential, power-dissipation, and relative-intensity utilities). For those, see the sections below and the demo notebooks.

## Getting Started

pyPI requires **Python 3.11+** to run.
To get pyPI up and running on your system, clone the repository and ensure that you have the required dependencies.

### Installation

[![PyPI version](https://badge.fury.io/py/tcpypi.svg)](https://badge.fury.io/py/tcpypi)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/tcpypi.svg)](https://anaconda.org/conda-forge/tcpypi)

Install from [PyPI](https://pypi.org/project/tcpypi/):

```bash
pip install tcpypi
```

or from [conda-forge](https://anaconda.org/conda-forge/tcpypi):

```bash
conda install -c conda-forge tcpypi
```

### tcpyPI Dependencies

Core (installed automatically): [NumPy](https://numpy.org/) (`>=1.26`) and [Numba](http://numba.pydata.org/) (`>=0.59`).

Optional extras:

```bash
pip install "tcpypi[xarray]"   # xarray + h5netcdf, to apply PI over gridded datasets
pip install "tcpypi[demo]"     # the above + matplotlib + pooch, for the example notebooks
```

Applying PI over large datasets is easy and highly recommended via [xarray](https://docs.xarray.dev/) (see the Quickstart below). tcpyPI's Numba kernels use on-disk caching (`cache=True`), so after the first call the compiled code is reused across sessions.

### Quickstart (5 minutes)

Compute PI for a single sounding (ordered surface → top, i.e. decreasing pressure):

```python
import numpy as np
from tcpyPI import pi

P  = np.array([1000.,975,950,925,900,850,800,750,700,600,500,400,300,250,200,150,100,70,50])  # hPa
TC = np.array([28.,25,24,23,22,19,16,13,11,5,-2,-11,-27,-37,-49,-65,-79,-73,-64])              # deg C
R  = np.array([18.,18,16,13,12,10,9,7,4,1.7,1.7,0.2,0.1,0.11,0.05,0.014,0.003,0.002,0.002])    # g/kg (mixing ratio)
SSTC, MSL = 30.0, 1010.0   # deg C, hPa

vmax, pmin, ifl, t0, otl = pi(SSTC, MSL, P, TC, R)
print(f"PI = {vmax:.1f} m/s, Pmin = {pmin:.1f} hPa, flag = {ifl}")
# PI = 81.9 m/s, Pmin = 903.3 hPa, flag = 1
```

Apply it over a gridded dataset with xarray (vectorized; this is what `run_sample.py` does):

```python
import xarray as xr
from tcpyPI import pi

ds = xr.open_dataset("data/sample_data.nc")   # dims: month, p, lat, lon
vmax, pmin, ifl, t0, otl = xr.apply_ufunc(
    pi, ds.sst, ds.msl, ds.p, ds.t, ds.r,
    kwargs=dict(CKCD=0.9),
    input_core_dims=[[], [], ["p"], ["p"], ["p"]],
    output_core_dims=[[], [], [], [], []],
    vectorize=True,
)
```

For very large grids, chunk the dataset with Dask and add `dask="parallelized", output_dtypes=[float]*5` to `apply_ufunc` — PI is computed independently per column, so it parallelizes trivially.

### Python Implementation of "pc_min" (BE02 PI Calculator)

[pi.py](src/tcpyPI/pi.py) is the Python function which directly computes PI given atmospheric and ocean state variables (akin to the BE02 algorithm MATLAB implementation [pc_min.m](matlab_scripts/pc_min.m)). Given input vector columns of environmental atmospheric temperatures (T) and mixing ratios (R) on a pressure grid (P), sea surface temperatures (SST), and mean sea-level pressures (MSL), the algorithm outputs potential intensity, the outflow level, the outflow temperature, and the minimum central pressure, and a flag that shows the status of the completed PI calculation. pyPI is an improvement on pcmin in that it handles missing values depending on user input flags.

Users who want to apply the PI calculation to a set of local environmental conditions need only to download [pi.py](./src/tcpyPI/pi.py), organize their data appropriately, and call the function to return outputs, e.g.:
```
(VMAX,PMIN,IFL,TO,LNB)=pi(SST,MSL,P,T,R)
```

### Sensitivity & Configuration Options

`tcpyPI.pi()` exposes a few key “sensitivity knobs” that can materially change results. The most commonly explored options are:

- `CKCD` (default `0.9`): exchange coefficient ratio (`Ck/Cd`).
- `diss_flag` (default `1`): include (`1`) or exclude (`0`) dissipative heating.
- `ascent_flag` (default `0`): reversible (`0`) vs pseudo-adiabatic (`1`) ascent assumption.
- `V_reduc` (default `0.8`): reduction from gradient wind to 10 m wind speed.
- `ptop` (default `50` hPa): top pressure bound; setting too high can bias outflow level/temperature.
- `miss_handle` (default `1`): missing-profile handling; `1` returns missing on any NaNs.
- `outflow_source` (default `"cape_star"`): how outflow level/temperature are defined:
  - `"cape_star"`: outflow is the LNB (`OTL`) and temperature (`T0`) from the saturated CAPE* calculation (default; BE02/pcmin behavior).
  - `"cape_env"`: outflow is the LNB and temperature from the environmental CAPE (CAPEenv) on the final convergence iteration (see Gilford et al. 2021 discussion).

Example:

```python
from tcpyPI import pi

vmax, pmin, ifl, t0, otl = pi(SSTC, MSL, P, TC, R, CKCD=0.9, diss_flag=1, outflow_source="cape_env")
```

### Log Decomposition (Wing et al. 2015)

tcpyPI provides a simple API for log-decomposing PI into efficiency, disequilibrium, and `Ck/Cd` terms (where `lnpi = ln(V^2)`):

```python
from tcpyPI import pi, pi_log_decomposition

vmax, pmin, ifl, t0, otl = pi(SSTC, MSL, P, TC, R, CKCD=0.9)   # SSTC/TC in Celsius
lnpi, lneff, lndiseq, lnCKCD = pi_log_decomposition(vmax, SSTC, t0, CKCD=0.9, sst_units="C")
```

`pi_log_decomposition` mirrors `vi_log_decomposition` and `gpi_log_decomposition`: give it the (already-computed) index inputs and it returns the additive log-space terms.

### Ventilation Index (Tang & Emanuel 2012)

tcpyPI includes an implementation of the Tang & Emanuel (2012) ventilation index:

Equations (TE12):

$$\Lambda = \left(\frac{V_{\mathrm{shear}}}{V_{PI}}\right)\,\chi_m$$

$$\chi_m = \frac{s_m^* - s_m}{s_{SST}^* - s_b}$$

where $V_{\mathrm{shear}}$ is the 850–200 hPa shear magnitude (m/s), $V_{PI}$ is potential intensity (m/s),
and $s_m^{\ast}, s_m, s_{SST}^{\ast}, s_b$ are moist entropies (e.g., $J\,kg^{-1}\,K^{-1}$) defined in TE12.

```python
import tcpyPI

# χm from a thermodynamic profile (TE12 Eq. 2)
chi_m = tcpyPI.entropy_deficit_te12_from_profile(
    P, TC, R, T_units="C", q_units="g/kg",
    SST=SSTC, SST_units="C", psfc_hpa=MSL,
    T2m=TC[0], q2m=R[0],
    entropy_method="emanuel94",                # or "bryan2008" (TE12-consistent)
    s_m_star_source="env_T_at_pmid",           # or "moist_adiabat_from_sst"
    sb_source="t2m",                           # or "layer_1000_900", "lowest_level"
    p_mid_hpa=600.0,
)

# Λ from shear and PI (TE12 Eq. 1)
Lambda = tcpyPI.ventilation_index(u_shear, VMAX, chi_m, formulation="te12")
```

Key options:
- `entropy_method`: `"emanuel94"` (default; consistent with most of tcpyPI) or `"bryan2008"` (matches TE12 Eq. 3).
- `s_m_star_source`: `"env_T_at_pmid"` (default) or `"moist_adiabat_from_sst"` (saturated parcel lifted from SST to `p_mid_hpa`).
- `sb_source`: `"t2m"` (default; uses `T2m/q2m` at `psfc_hpa`), with fallbacks available via `sb_fallback`.

See `notebooks/ventilation_index_demo.ipynb` for a step-by-step worked example and visualization.

### Power Dissipation Index (PDI)

tcpyPI includes a Power Dissipation Index (PDI) utility for storm-track time series:

Equations (Emanuel 2005):

$$\mathrm{PDI} = \sum_{t} V_{\max}(t)^3\,\Delta t$$

where $V_{\max}$ is the maximum sustained wind speed and $\Delta t$ is the time step. The
`"e05_si"` formulation returns the SI-consistent sum (e.g., $m^3\,s^{-2}$ if $V$ is m/s and $\Delta t$ is s),
and `"e05_1e11"` returns `PDI / 1e11` for plotting convenience.

```python
import tcpyPI

# Vmax time series (example: m/s) and time step (example: hours)
PDI = tcpyPI.power_dissipation_index(
    vmax_series,
    dt_hours,
    wind_units="m/s",    # or "kt"
    dt_units="h",        # or "s"
    formulation="e05_si" # or "e05_1e11" (scaled for plotting)
)
```

Options:
- `nan_policy`: `"propagate"` (default) returns NaN if any timestep is missing; `"omit"` drops missing contributions.
- `dt` may be a scalar (constant timestep) or an array broadcastable to `vmax_series` (variable timestep).

See `notebooks/power_dissipation_index_demo.ipynb` for an end-to-end example using observed intensity from IBTrACS (Hurricane Milton, 2024) and computing PDI.

### Relative Intensity (ν)

tcpyPI provides a small helper for **relative intensity**, commonly defined as:

$$\nu(t) = V_{\max}(t) / V_{PI}(t)$$

where $\nu$ is unitless, $V_{\max}$ is the storm intensity, and $V_{PI}$ is the local potential intensity.

Example:

```python
import tcpyPI

nu = tcpyPI.relative_intensity(vmax_series, vpi_series)
```

### Genesis Potential Index (GPI)

tcpyPI includes a Genesis Potential Index (GPI) utility:

Equations (as implemented here):

For `"en04"` (Emanuel & Nolan 2004; Camargo et al. 2007):

$$\mathrm{GPI} = \left(10^5\,|\eta|\right)^{3/2}\,\left(\frac{\mathrm{RH}}{50}\right)^3\,\left(\frac{V_{PI}}{70}\right)^3\,\left(1 + 0.1\,V_{\mathrm{shear}}\right)^{-2}$$

For `"e10"` (Emanuel 2010, JAMES, Eq. 11):

$$\mathrm{GPI} = |\eta|^{3}\,\chi_m^{-4/3}\,\max(V_{PI}-35,\,0)^{2}\,\left(25 + V_{\mathrm{shear}}\right)^{-4}$$

where $|\eta|$ is low-level absolute vorticity ($s^{-1}$; commonly 850 hPa) and $V_{PI}$ is potential intensity (m/s). For `en04`, RH is midlevel relative humidity (%; commonly 700 hPa) and $V_{\mathrm{shear}}$ is the 850–200 hPa shear magnitude (m/s). For `e10`, $\chi_m$ is the nondimensional midlevel (600 hPa) moist-entropy deficit of Tang & Emanuel (2012) (i.e. `tcpyPI.entropy_deficit_te12_from_profile`), and $V_{\mathrm{shear}}$ is the 850–250 hPa shear magnitude (m/s).

**GPI is a relative index** (defined up to a constant of proportionality), so `en04` and `e10` magnitudes are *not* directly comparable to each other.

```python
import tcpyPI

# Inputs (typical choices):
# - abs_vort: absolute vorticity at ~850 hPa (s^-1); the sign is ignored (|eta| is used)
# - v_shear: deep-layer shear magnitude (m/s)
# - v_pot: potential intensity (m/s), e.g. from tcpyPI.pi()

# en04 (RH-based); rh_mid = RH at ~700 hPa (%)
gpi_en04 = tcpyPI.genesis_potential_index(abs_vort, rh_mid, v_shear, v_pot, formulation="en04")

# e10 (entropy-deficit-based): compute chi_m via the ventilation-index utility
# (see the Ventilation Index section for the full call and options), then:
gpi_e10 = tcpyPI.genesis_potential_index(
    abs_vort, v_shear=v_shear, v_pot=v_pot, formulation="e10", chi=chi_m
)
```

Formulations:
- `"en04"`: Emanuel & Nolan (2004) / Camargo et al. (2007), relative-humidity-based
- `"e10"`: Emanuel (2010), midlevel moist-entropy-deficit-based ($\chi_m$); requires `chi`

See `notebooks/genesis_potential_index_demo.ipynb` for an end-to-end worked example.

### Running a pyPI Sample

Included in the pyPI release is a sample script [run_sample.py](run_sample.py) which runs a compact sample dataset (ERA5 October 2024 monthly-mean conditions; “Milton era”) through `pi.py`, vectorizes the output, and performs several simple analyses. To run, simply:
```
python run_sample.py
```
and examine the outputs locally produced in [full_sample_output.nc](./data/full_sample_output.nc).

### Use tcpyPI with AI (MCP server)

tcpyPI ships a [Model Context Protocol](https://modelcontextprotocol.io/) server that gives AI assistants (Claude Desktop, Claude Code, Cursor, and other MCP clients) direct access to the validated PI calculation and genesis diagnostics — instead of a model improvising moist thermodynamics from memory. Install and register:

```
pip install "tcpypi[mcp]"
```

- **Claude Code:** `claude mcp add tcpypi -- tcpypi-mcp`
- **Claude Desktop / Cursor / other clients:** add a stdio server entry running `tcpypi-mcp` (or, without a prior install, command `uvx` with args `--from "tcpypi[mcp]" tcpypi-mcp`).

Six tools are exposed: `compute_pi` (single sounding), `compute_pi_grid` (netCDF file in → netCDF file out, whole-field fast path), `decompose_pi` (Wing et al. 2015 log decomposition), `ventilation_index` (TE12), `genesis_potential_index` (EN04/E10), and `power_dissipation_index` (Emanuel 2005, from intensity time series). Every result carries provenance (package version, DOI, citations, knob settings) and the `ifl` status contract (only `ifl == 1` output is trustworthy). Units are fixed (°C, hPa, g/kg) and **never converted silently**: gridded inputs must carry matching netCDF `units` attributes, and mismatches are rejected with a clear error.

## File Descriptions

#### Key files
* **[pi.py](./src/tcpyPI/pi.py)** - The primary function of pyPI, that computes and outputs PI (and associated variables) given atmospheric and ocean state variables.
* **[run_sample.py](run_sample.py)** - Example script that computes PI and accompanying analyses over the entire sample dataset

#### Data
* [sample_data.nc](./data/sample_data.nc) - Sample state variables from ERA5 2024 monthly means over the North Atlantic (12 months; ~2° subsample; SST masked to ocean). Includes the core PI inputs (SST, MSL, T, and the true water-vapor mixing ratio `r = q/(1−q)` converted from ERA5 specific humidity) plus winds, RH, and vorticity for the GPI/ventilation-index analyses. Built by [data/build_sample_data_era5_oct2024.py](./data/build_sample_data_era5_oct2024.py) from the (untracked, regenerable) ERA5 downloads.
* [era5_demo_subset.nc](./data/era5_demo_subset.nc) - Small single-column ERA5 environment (Hurricane Milton 2024 reference location) used by the ventilation-index and GPI demo notebooks.
* [mdr.json](./data/mdr.json) - Main Development Region definitions from [Gilford et al. (2017)](https://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-16-0827.1)
* [raw_sample_output.nc](./data/raw_sample_output.nc) - Sample outputs from pi.py *only*, created by run_sample.py
* [full_sample_output.nc](./data/full_sample_output.nc) - Full set of sample outputs from pi.py plus analyses (efficiency, disequilibrium, log decomposition)
* [matlab_pi_reference_2024.nc](./data/matlab_pi_reference_2024.nc) - BE02 MATLAB (pc_min) reference outputs for the 2024 sample, generated by [reference_calculations.m](./matlab_scripts/reference_calculations.m) and used by verify_pi.ipynb

#### Validation and Testing Notebooks
* **[test_pi_calc.ipynb](./notebooks/test_pi_calc.ipynb)** - Single call of pi.py on the 2024 sample plus a quick run-time check
* **[verify_pi.ipynb](./notebooks/verify_pi.ipynb)** - Validates pyPI outputs against the BE02 MATLAB algorithm (pc_min); reference outputs for the 2024 sample are generated by [reference_calculations.m](./matlab_scripts/reference_calculations.m)
* **[illustrate_numerical_instability.ipynb](./notebooks/illustrate_numerical_instability.ipynb)** - Documents a marginal (near-neutral) sounding where the BE02 pressure solver does not converge (`IFL=2`), contrasted with an adjacent-day profile that does; explains the fixed-point oscillation and the `IFL`-filtering guidance
* **[sample_output_analyses.ipynb](./notebooks/sample_output_analyses.ipynb)** - North Atlantic 2024 seasonal analysis: PI, efficiency/disequilibrium, GPI, and ventilation index, each with a log-decomposition
* **[efficiency_demo.ipynb](./notebooks/efficiency_demo.ipynb)** - The tropical-cyclone (Carnot) efficiency term $(T_s-T_0)/T_0$ from BE02: what it represents and its sensitivity to SST and outflow temperature
* **[dissipative_heating_effect.ipynb](./notebooks/dissipative_heating_effect.ipynb)** - Compares PI with and without dissipative heating (`diss_flag`); reproduces the ~20–25% intensification noted by Bister & Emanuel (1998)
* **[ventilation_index_demo.ipynb](./notebooks/ventilation_index_demo.ipynb)** - TE12 ventilation index: end-to-end calculation and diagnostics
* **[power_dissipation_index_demo.ipynb](./notebooks/power_dissipation_index_demo.ipynb)** - PDI: computed from observed IBTrACS intensity for Hurricane Milton (2024)
* **[genesis_potential_index_demo.ipynb](./notebooks/genesis_potential_index_demo.ipynb)** - GPI (`en04` and Emanuel-2010 `e10`): worked example and sensitivities

#### Misc.
* **[utilities.py](./src/tcpyPI/utilities.py)** - Set of functions used in the pyPI codebase
* **[constants.py](./src/tcpyPI/constants.py)** - Set of meteorological constants used in the pyPI codebase
* **[vi.py](./src/tcpyPI/vi.py)** - TE12 ventilation index and entropy deficit utilities
* **[pdi.py](./src/tcpyPI/pdi.py)** - Power dissipation index utility
* **[gpi.py](./src/tcpyPI/gpi.py)** - Genesis potential index utility
* **[reference_calculations.m](./matlab_scripts/reference_calculations.m)** - Reads `data/sample_data.nc` and generates BE02 MATLAB (pc_min) reference outputs used by verify_pi.ipynb
* **[pc_min.m](./matlab_scripts/pc_min.m)** - Original BE02 algorithm from MATLAB, adapted and used to produce analyses of Gilford et al. ([2017](https://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-16-0827.1); [2019](https://journals.ametsoc.org/doi/10.1175/MWR-D-19-0021.1))
* **[clock_pypi.ipynb](./notebooks/clock_pypi.ipynb)** - Notebook estimating the time it takes to run pyPI on a laptop

## Author

* **Daniel M. Gilford, PhD** - *Creation, Development, & Maintenance* - [GitHub](https://github.com/dgilford)

### Contributor(s)
* **Ben Mares** - *Modernization* - [GitHub](https://github.com/maresb)
* **Daniel Rothenberg, PhD** - *Numba Optimization & Sample Code* - [GitHub](https://github.com/darothen)

## Development & AI transparency

tcpyPI v2.0 was developed in part with Claude Code and carefully vetted by me. Any errors or oversights are my own.

Recent maintenance of tcpyPI (v1.4.x) has been **AI-assisted** (with Claude Code) for code review, refactoring, documentation, and test scaffolding. Every change is **human-reviewed and hand-edited by the author** and is guarded by the automated test suite: the `pi()` doctests pin scalar outputs to ~13 significant figures, and `tests/test_run_sample.py` pins the gridded sample outputs at `rtol=1e-13`. AI assistance is used as a productivity tool — the science, defaults, and interfaces remain the author's.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Kerry Emanuel (MIT) - Development of potential intensity theory; encouragement and permission to pursue Python implementation
* Susan Solomon (MIT), Paul O'Gorman (MIT), Allison Wing (FSU) - Helpful conversations, advice, and suggestions on TC PI research
* Dan Chavas (Purdue), Jonathan Lin (MIT), Raphael Rousseau-Rizzi (MIT) - Feedback on pyPI code and documentation
