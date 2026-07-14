# Repository Guidelines

Guidance for AI agents and human contributors working in this repo. (`CLAUDE.md` points here.)

## Project Structure & Module Organization
Core library in `src/tcpyPI/`:
- `pi.py` — the BE02 PI algorithm (`pi()`, `cape()`, entropy solver), the whole-field fast path `pi_field()` (parallel guvectorize; ~2 orders of magnitude over per-column loops), and the Wing et al. (2015) log-decomposition `pi_log_decomposition()`. Numba-compiled. `cape()` keeps the original pcmin.m per-call structure (an environment-setup hoist was certified bit-identical but reverted: negligible speed gain did not justify structural drift from the MATLAB reference), with one deliberate algorithmic divergence: CAPE is the *maximum* running buoyancy work over candidate outflow levels (max-W), not pcmin.m's integrate-to-the-highest-positive-level rule, whose discontinuity caused the issue #77 solver failures (see `discontinuity_analysis/`).
- `utilities.py` — thermodynamic helpers, unit/longitude conversions, PI efficiency/disequilibrium (`decompose_pi` scalar core).
- `constants.py` — meteorological constants (fully documented, incl. the deliberate modified `CL=2500`); `numba.py` — the `@njit`/`guvectorize` wrapper (graceful pure-Python fallback when Numba is absent).
- `vi.py` — TE12 ventilation index + entropy-deficit (`chi_m`) + `vi_log_decomposition`. `gpi.py` — genesis potential index (`en04`, `e10`) + `gpi_log_decomposition`. `pdi.py` — power dissipation index. `relative_intensity.py` — ν = V/V_PI. (VI/GPI/PDI are experimental.)

Sample workflow in `run_sample.py` (writes `data/raw_sample_output.nc`, `data/full_sample_output.nc`). Tests in `tests/`. Notebooks in `notebooks/`; figures in `figures/`; legacy BE02 MATLAB in `matlab_scripts/`.

## Build, Test, and Development Commands
- `pip install -e .` (Python 3.11+). Optional extras: `.[xarray]`, `.[demo]`.
- Tests and tooling use pixi. **Note:** a stale `pixi` shim may shadow the real binary — use `~/.pixi/bin/pixi` if `pixi` fails with a module error.
  - `pixi run -e test-latest basic-tests` — core suite + doctests.
  - `pixi run -e test-xarray-latest xarray-tests` — sample-output regression (rtol=1e-13 pins).
  - `pixi run -e test-min basic-tests` — minimum supported floor (py3.11 / numpy 1.26 / numba 0.59).
- Lint/format/type-check (enforced by the CI `lint` job on every push/PR):
  - `pixi run -e lint lint` (ruff check) · `pixi run -e lint format-check` / `format` (ruff format) · `pixi run -e lint typecheck` (mypy).
  - Contributors can also `pip install pre-commit && pre-commit install` (ruff + ruff-format hooks).
- Data: `pixi run download-era5` (needs a CDS API key in `~/.cdsapirc`) then `pixi run build-sample` regenerates `data/sample_data.nc` + `data/era5_demo_subset.nc`. `python run_sample.py` regenerates the output artifacts.
- Benchmarks: the `bench` pixi env adds Skyborn for the optional comparison in `notebooks/clock_pypi.ipynb`. To run that notebook in JupyterHub, register the env as a kernel: `pixi run -e bench python -m ipykernel install --user --name tcpypi-bench --display-name "tcpyPI (bench)"`.

## Data Policy
- Track only the small samples (`sample_data.nc` ~6 MB, `era5_demo_subset.nc`, run_sample outputs, the MATLAB reference). `.gitignore` uses an explicit allowlist — do **not** commit raw ERA5 downloads or other large `.nc`.
- The tracked sample is regional (North Atlantic 2024), subsampled to ~2° (decimation, no interpolation), SST masked to ocean; humidity is stored as the **true water-vapor mixing ratio `r`** (converted from ERA5 specific humidity via `r = q/(1−q)`), with ERA5 relative humidity as `rh`. Regenerate a wider/global sample by widening the `area` in `data/download_era5.py`.

## Coding Style & Scientific Invariants
- PEP 8 via ruff (line length 100); snake_case. Document units for scientific quantities; prefer keyword-only args for flags on new APIs (see `pi_field`).
- **Never enable Numba `fastmath`** anywhere: the PI/CAPE missing-value logic relies on IEEE NaN semantics (`NaN > 0` is False, min/max NaN propagation). `@njit` already defaults to `cache=True` (see `numba.py`); don't remove it.
- **Kernel-change policy (bit-for-bit identicality):** the compiled PI/CAPE kernels implement the published GMD-2021/pcmin algorithm — plus one documented, deliberate divergence (max-W CAPE, issue #77; validated in `discontinuity_analysis/`: no flag changes and 98.1% of VMAX within 0.1 m/s on the 2024 sample) — and any *other* change to them must be *math-by-math identical* — proven, not asserted. Acceptable transforms: selections (scalar `min`/`max` with preserved argument order), hoisting/deduplicating pure computations, memory-layout changes. **Never** reorder/reassociate arithmetic. Certify with the byte-equality harness (`.ai/bit_identity_harness.py`: old-commit worktree vs working tree over the full sample × option combos through both entry points, direct `cape()` calls, and a 50k fuzz corpus), and re-certify on the final post-`ruff format` bytes.
- `pi()`/`pi_field()` expect SST/T in °C, MSL/P in hPa, mixing ratio in g/kg; any vertical profile order is accepted (sorted internally to decreasing pressure). NaN P/MSL → `IFL=3`; NaN SST → `IFL=0`; period-2 solver oscillations with amplitude ≲1 hPa are rescued to the cycle midpoint (`IFL=1`; documented divergence from pcmin.m, which returns missing for those soundings). Flags tripped by the environmental-CAPE call persist (#78), so `IFL∈{0,3}` can accompany finite outputs — **`IFL==1` is the trustworthy-output gate**.
- Keep changes test-pinned: `pi()` doctests pin outputs to ~13 sig figs; `tests/test_run_sample.py` pins the sample outputs at `rtol=1e-13`; `tests/test_oscillation_rescue.py` pins bit-identity of a converging control profile; `tests/test_pi_field.py` pins `pi_field == pi` exact equality.

## Commit & Pull Request Guidelines
- Imperative commit summaries; group related changes. Reference issues/PRs; note any data or API changes.
- Before pushing, run `basic-tests`, the lint tasks, and (when touching `run_sample.py`/data) `xarray-tests`. Do not commit large binaries.
