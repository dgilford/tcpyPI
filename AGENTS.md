# Repository Guidelines

Guidance for AI agents and human contributors working in this repo. (`CLAUDE.md` points here.)

## Project Structure & Module Organization
Core library in `src/tcpyPI/`:
- `pi.py` — the BE02 PI algorithm (`pi()`, `cape()`, entropy solver). Numba-compiled.
- `utilities.py` — thermodynamic helpers, unit/longitude conversions, PI efficiency/disequilibrium.
- `constants.py` — meteorological constants; `numba.py` — the `@njit` wrapper.
- `analyses.py` — Wing et al. (2015) PI log-decomposition API.
- `vi.py` — TE12 ventilation index + entropy-deficit (`chi_m`). `gpi.py` — genesis potential index (`en04`, `e10`). `pdi.py` — power dissipation index. `relative_intensity.py` — ν = V/V_PI. (VI/GPI/PDI are experimental.)

Sample workflow in `run_sample.py` (writes `data/raw_sample_output.nc`, `data/full_sample_output.nc`). Tests in `tests/`. Notebooks in `notebooks/`; figures in `figures/`; legacy BE02 MATLAB in `matlab_scripts/`.

## Build, Test, and Development Commands
- `pip install -e .` (Python 3.10+). Optional extras: `.[xarray]`, `.[demo]`.
- Tests use pixi. **Note:** a stale `pixi` shim may shadow the real binary — use `~/.pixi/bin/pixi` if `pixi` fails with a module error.
  - `pixi run -e test-latest basic-tests` — core suite + doctests (49 tests).
  - `pixi run -e test-xarray-latest xarray-tests` — sample-output regression (2 tests).
  - `pixi run -e test-min basic-tests` — minimum supported floor (py3.10 / numpy 1.26 / numba 0.59).
- Data: `pixi run download-era5` (needs a CDS API key in `~/.cdsapirc`) then `pixi run build-sample` regenerates `data/sample_data.nc` + `data/era5_demo_subset.nc`. `python run_sample.py` regenerates the output artifacts.

## Data Policy
- Track only the small samples (`sample_data.nc` ~6 MB, `era5_demo_subset.nc`, run_sample outputs). `.gitignore` uses an explicit allowlist — do **not** commit raw ERA5 downloads or other large `.nc`.
- The tracked sample is regional (North Atlantic), subsampled to ~2° (decimation, no interpolation), SST masked to ocean. Regenerate a wider/global sample by widening the `area` in `data/download_era5.py`.

## Coding Style & Scientific Invariants
- PEP 8, 4-space indent, snake_case. Document units for scientific quantities; prefer keyword args for flags.
- **Never enable Numba `fastmath`** anywhere: the PI/CAPE missing-value logic relies on IEEE NaN semantics (`NaN > 0` is False, min/max NaN propagation). `@njit` already defaults to `cache=True` (see `numba.py`); don't remove it.
- `pi()` expects SST/T in °C, MSL/P in hPa, mixing ratio in g/kg, profiles ordered surface→top (the wrapper now sorts to descending pressure, so any order is accepted).
- Keep changes test-pinned: `pi()` doctests pin outputs to ~13 sig figs; `tests/test_run_sample.py` pins the sample outputs at `rtol=1e-13`.

## Commit & Pull Request Guidelines
- Imperative commit summaries; group related changes. Reference issues/PRs; note any data or API changes.
- Run `basic-tests` (and `xarray-tests` when touching `run_sample.py`/data) before pushing. Do not commit large binaries.
