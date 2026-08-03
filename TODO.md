# tcpyPI to-do

Tracked work items for the v2.0 release. Durable "how this repo works" knowledge belongs in
[AGENTS.md](AGENTS.md); this file is for things that still need doing.

## v2.0

### Reconcile longitude conventions between `sample_data.nc` and `mdr.json`

The two tracked data files disagree on longitude convention:

- `data/sample_data.nc` is on **0–360** (250 → 320 for the North Atlantic sample)
- `data/mdr.json` basin boxes are on **−180/180** (e.g. `na` = −95 → −50)

Selecting an MDR box straight out of the JSON against the sample returns an **empty selection
with no error** (`lon` dim of size 0), silently producing all-NaN results instead of failing.
Currently documented as a caveat in AGENTS.md (Data Policy); v2.0 should actually fix it.

Options, with their costs:

1. **Convert `mdr.json` to 0–360** to match the data. Cheapest. Check first whether anything else
   consumes `mdr.json` and whether changing it muddies the Gilford et al. (2017) region-definition
   provenance.
2. **Rebuild `sample_data.nc` on −180/180.** Most invasive — requires regenerating the run_sample
   outputs and the MATLAB reference, which touches the `rtol=1e-13` regression pins.
3. **Add wrap-aware region selection** to tcpyPI so both conventions work. Check `~/doombot` for an
   existing longitude-wrapping helper before writing a new one.

Whichever is chosen, **add a test asserting that an `mdr.json` box selected against
`sample_data.nc` returns a non-empty subset**, so the silent-empty-selection failure mode cannot
regress.

*Found 2026-08-03 while building the AMS 2027 abstract figure.*
