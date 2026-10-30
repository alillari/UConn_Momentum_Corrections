# Reusable exclusivity-window workflow

The workflow derives **candidate** exclusivity windows in explicit
sector-separated 3D cells. It does not derive a momentum correction. The fitting
library accepts an ordinary 1D count histogram and knows nothing about ROOT,
runs, particles, sectors, filenames, or histogram-title conventions.

## Components and boundaries

| Module | Responsibility |
| --- | --- |
| `models.py` | Bin-integrated signal/background components and analytic derivatives |
| `peakfit.py` | Poisson fitting, quality checks, model/range comparison |
| `cuts.py` | Serializable windows and vectorized decisions for supplied cell IDs |
| `binning.py` | Coordinate convention, exact 3D boundaries and cell IDs |
| `kinematics.py`, `skim.py` | Physics adapter shared by diagnostics and cut development |
| `root_skim.py` | Logical branch mappings and stable ROOT row access |
| `calibration.py` | Cache construction, run separation and batch fitting |
| `evaluation.py` | Applying windows, friend masks and retention accounting |
| `cut_report.py` | Rendering saved results; no fit or cut policy |
| `recipe.py` | Strict configuration validation |
| `channels.py` | Charged identities, measured topologies and four-vector observables |

No model or cut is chosen by a hard-coded run, sector, cell ID, or plot title.
Numerical cut policy, fit families/ranges, admission counts, stability gates,
and run partitions are in the recipe. New signal/background components extend
the component registry; ROOT iteration and cut application do not change.
The current signal interface has a center and width, with bin probabilities
and derivatives; bundled backgrounds are positive Bernstein polynomials.

The schema-2 workflow supports electron, pi_plus, pi_minus and proton tracks.
It includes epipn, fully measured elastic ep, and fully measured ep pi+ pi-
templates. Negative pions require all four tracks measured. See
[channel development](channels.md) for the input/observable contract and the
present FD-only boundary. Original schema-1 RG-K commands and artifacts remain
supported. No measured elastic or double-pion dataset has been validated here.

## Commands

Install `tools/requirements-cuts.txt` into the chosen Python environment. The
full example uses the RG-K recipe:

```bash
python tools/exclusivity_cuts.py run \
  --input /path/to/prepare_root_files \
  --diagnostics-config CLAS12/JSON_configs/rgk_epip_diagnostics.json \
  --binning-config CLAS12/JSON_configs/rgk_epip_binning.json \
  --recipe CLAS12/JSON_configs/rgk_epip_cuts.json \
  --output /path/to/exclusivity_output \
  --workers 4
```

Stages can be run independently:

```bash
python tools/exclusivity_cuts.py cache --input /path/to/data \
  --diagnostics-config CLAS12/JSON_configs/rgk_epip_diagnostics.json \
  --binning-config CLAS12/JSON_configs/rgk_epip_binning.json \
  --recipe CLAS12/JSON_configs/rgk_epip_cuts.json --output /path/to/cache

python tools/exclusivity_cuts.py derive --cache /path/to/cache \
  --recipe CLAS12/JSON_configs/rgk_epip_cuts.json \
  --output /path/to/model --workers 4

python tools/exclusivity_cuts.py apply --input /path/to/data \
  --selection /path/to/model/cuts.json --output /path/to/evaluation

python tools/exclusivity_cuts.py report --cache /path/to/cache \
  --model /path/to/model --evaluation /path/to/evaluation \
  --output /path/to/report
```

`derive` never reads validation spectra. Changing model or window policy does
not require rereading the skim; changing partitions, roles or binning does
require a corresponding new cache. The tests explicitly verify that changing
the cached validation spectrum cannot change fitted windows. Cache and
selection application are usable without SciPy; optimization and reporting
require it. `--no-masks` skips friend-file export for exploratory evaluations.

## Fitting and cut policy

Fits use an unweighted Poisson count likelihood, including empty bins. The
expected signal and background are integrated over actual bin edges. The
Gaussian signal yield is defined over its full PDF; each positive background
basis component is normalized over the fit domain. Analytic derivatives are
checked against finite differences, and observable-unit scaling is tested.
Poisson histogram fitting is also described in the [ROOT documentation](https://root.cern.ch/doc/master/classTH1.html).

For each configured fit range, the fitter tries the model family and width
initializations. It checks convergence, local covariance, residual deviance,
peak/background likelihood gain, proximity to parameter bounds and parameter
uncertainties. The likelihood gain is a local algorithmic contrast diagnostic,
not a discovery significance after scanning models/centers.

AICc is compared **only within the same fit domain**. Plausible models within
the configured AICc difference and winners from different domains must pass
center/width stability gates. One actual fit supplies the candidate window.
Correlated alternatives are never averaged to shrink statistical errors. Their
spread is a sensitivity diagnostic, not a calibrated systematic uncertainty.
The local Fisher covariance is likewise conditional on the chosen model.

The initial window is the fitted center minus/plus three fitted widths. The
lower/upper multiples are separately configurable. A window extending beyond
the fitted domain is marked for review rather than extrapolated or clipped.
No hand-adjusted special cases were added for individual data cells.

## Artifacts and application semantics

- `cache/cache.json`, `manifest.json`, `histograms.npz`: source hashes,
  configuration snapshots, run-separated spectra, and complete accounting.
- `model/cuts.json`: every cell, bounds, candidate window, status, reasons and
  provenance. `ready` means only that configured algorithmic gates passed.
- `model/fits.json.gz`: every attempted model/range, its parameters, quality
  diagnostics and rejection reasons. `*_windows.csv` provides a review table.
- `evaluation/masks/*.root`: one ROOT friend per source file, in unchanged
  entry order, with run/event IDs and separate per-role cell/pass/reason fields.
- `evaluation/evaluation.json`, `evaluation_histograms.npz`: support losses,
  cut rejection, run comparisons, width sensitivity, display losses and hashes.

The window table is separate from coordinate transforms and observables:

```python
from momcorr.cuts import WindowTable

table = WindowTable(selection["cells"][role], grid.size)
passed, reason = table.apply(cell_ids, missing_mass)
```

Each particle role is evaluated independently; there is no implicit AND/OR.
Decision codes are 0 pass, 1 outside grid, 2 nonfinite observable, 3 unsupported
cell, and 4 outside the window. Unsupported cells have no borrowed or default
window. A false pass flag with reason 3 is an absence of calibration coverage,
not a measured physics rejection. Input files and their momenta are preserved.

## Validation limits

Before/after missing mass is a selection diagnostic on the selected variable.
The inferred pion momentum residual is correlated with missing mass, so its
narrowing does not independently establish correction performance. Pass
fractions are retention, not signal efficiency or purity. The configured
run-held-out evaluation is exploratory: binning was explored on the full sample
before the fit split. Further independent physics validation is required before
treating the candidate as a production selection.

Unit/integration tests cover generated physical closure, model derivatives,
injected peak recovery, a background-only spectrum, observable-unit scaling,
empty bins, exact 3D edges, run separation, validation isolation, unsupported
cells, friend-mask alignment and unchanged inputs. A separate 80-toy check of
the correctly specified Gaussian/sloped-background model gave 67.5% coverage
for nominal 68% center intervals; it does not validate arbitrary backgrounds.
