# Momentum Corrections

This repository contains ongoing development for momentum and energy correction studies for CLAS12 data.  
⚠️ **Work in progress** — this code is under active development and is **not intended for production use yet**.

## Purpose

The goal of this project is to provide a modular and configurable framework for:
- Applying energy and momentum corrections to reconstructed particles
- Computing derived quantities
- Producing diagnostic histograms using ROOTDataFrame

The design emphasizes **reusability**, **clear separation of responsibilities**, and **dataset/particle-level configurability**.

## Project Structure

- `include/`  
  Core headers and classes (e.g. particle abstraction, correction interfaces, configuration structs).

- `libs/`  
  Implementation files for correction models and helper utilities.

- `CLAS12/analysis/`  
  Analysis executables and scripts using the core library.

- `CLAS12/JSON_configs/`  
  Temporary home for json files used as input to analysis files.

- `CMakeLists.txt`  
  Build configuration for compiling libraries and analysis code.

## Build

This project uses CMake and expects ROOT to already be available in the environment.

```bash
cd UConn_Momentum_Corrections
mkdir -p build
cd build
cmake ..
cmake --build .
```

## Executables

Executables are stored in CLAS12/executable/ .

## Initial Histogram Production

Histogram production now uses a strict v2 JSON config. The config defines:
- input tree name
- particle roles and branch names
- selected C++ observables
- histogram requests and detector binning
- output directory

To make initial histograms, run as such:
```bash
cd CLAS12/executable

./make_initial_histograms \
  /path/to/input/root/files \
  all \
  ../JSON_configs/pip_CD_config.json
```

changing the json file as needed.

The current example config, `CLAS12/JSON_configs/pip_CD_config.json`, declares
the `epip(N)` missing-mass and pion `deltaP` histogram workflow. Unknown schema
keys, missing input branches, invalid binning, unknown detector adapters, and
unknown observables are treated as fatal configuration errors.

## Current Framework Boundary

This branch establishes a durable histogram and observable foundation. It does
not yet specify how momentum corrections are derived or applied. Correction
production should plug into the provenance-rich histogram/observable outputs in
a later design pass once the correction method is defined.

## Uncut single-pion diagnostics (optional Python frontend)

`tools/exclusivity_diagnostics.py` reads two-track ROOT skims without building
the C++ executable. Install `tools/requirements-diagnostics.txt` in a Python
environment and run:

```bash
python tools/exclusivity_diagnostics.py \
  --input /path/to/prepare_root_files \
  --config CLAS12/JSON_configs/rgk_epip_diagnostics.json \
  --output /path/to/diagnostic_output
```

The example is specifically for RG-K fall 2018, 6535 MeV, Pass 2. Its nominal
6.53536 GeV beam energy and LH2 target have documented sources in the JSON.
The nominal energy is not a per-run RCDB measurement. The input has no
energy-loss corrections according to the data provider; other preparation
selections are not documented in the ten-branch skim.

This frontend applies no physics cuts, corrections, or candidate deduplication.
It produces `diagnostics.pdf`, PNG figures, `summary.json`, a SHA-256 input
`manifest.json`, and `histograms.npz`. The summary records numerical-invalid
rows and display-range losses, schema/run checks, repeated run/event pairs,
software versions, and code provenance. Repeated pairs are reported without
assuming they are duplicate events: multiple pion candidates can share an ID.

Missing mass squared retains its sign; the mass plot uses a signed square root.
The pion momentum constraint is diagnostic only, and ambiguous physical roots
are reported instead of choosing the root closest to the measured momentum.
Input sector labels are preserved without assuming a detector region. Missing
PID, vertex, timing, and track-quality quantities cannot be retuned from this
skim. The diagnostics JSON has its own schema and is not a C++ v2 config.

Run its independent physics and ROOT accounting checks with:

```bash
python -m unittest discover -s tests -p 'test_exclusivity_diagnostics.py'
```

To adjust report layout without rereading the data, use the same command with
`--render-only`. All numerical results remain in the saved summary/histograms.

## 3D binning before exclusivity fits

The binning frontend establishes **sector-separated cells in momentum magnitude,
theta, and shifted phi** for each particle. It compares coarse, intermediate
(`balanced`), and fine grids on the full input, and saves the missing-mass
spectrum of every cell in the selected grid. No fitting or selection is applied.

```bash
python tools/exclusivity_binning.py \
  --input /path/to/prepare_root_files \
  --diagnostics-config CLAS12/JSON_configs/rgk_epip_diagnostics.json \
  --binning-config CLAS12/JSON_configs/rgk_epip_binning.json \
  --output /path/to/3d_binning_output
```

The coordinate conventions, units, candidate grids, and selected schemes are
explicit in the binning JSON. The user-confirmed shifted-phi equations are
`phi_local - 30/p` for electrons and `phi_local + 32/(p - 0.05)` for positive
pions, with momentum in GeV/c and azimuth in degrees. `phi_local` uses Richard's
sector unwrap and subtracts `(sector - 1)*60`. The result is not wrapped again.
Theta is in degrees, independently of the older C++ branch conventions.

The default intermediate grids use 0.25 GeV/c momentum bins in the populated
region, variable-width theta bins, and 10-degree shifted-phi bins. Tails use
wider bins. The 0.85-1.05 GeV mass band is only a population diagnostic and
contains background; its counts do not establish signal yield or fit quality.
See `docs/rgk_3d_binning.md` for the measured tradeoffs and domain coverage.

Outputs are `binning_report.pdf`, PNG views, `binning.json`, `summary.json`,
`manifest.json`, `cells.npz`, and one cell-boundary/count CSV per particle.
`cells.npz` stores selected-grid mass spectra, explicit spectrum overflow
counts, every candidate grid's populations, and selected-grid run populations.
All sectors remain separate. Plot projections may sum over sectors, but the
saved cells never do. Empty and out-of-domain cells/rows are explicitly counted.

The Python cell-assignment API is in `python/momcorr/binning.py`. IDs are C-order
flat indices of `(sector, p, theta, shifted_phi)`. Intervals are half-open with
the final right edge included. Retrieve a cell's spectrum using the exact saved
config and `cells.npz`:

```python
import json
import numpy as np
from momcorr.binning import load_binning  # PYTHONPATH=/path/to/repo/python

config = json.load(open("binning.json"))
grids = load_binning(config)
role = "pion"
scheme = config["particles"][role]["selected_scheme"]
grid = grids[role][scheme]
cell_id = 100
with np.load("cells.npz") as cells:
    bounds = grid.bounds(cell_id)
    counts = cells[f"{role}_spectra"][cell_id]
    mass_edges = cells["mass_edges"]
```

Use `--render-only` to redraw an existing bundle; changed configurations require
a new analysis. Run all optional Python tests with `python -m unittest discover
-s tests`.

## Candidate exclusivity windows

The reusable statistical engine, explicit cut table, and run-separated ROOT
workflow are described in [docs/exclusivity_workflow.md](docs/exclusivity_workflow.md).
Install `tools/requirements-cuts.txt`, then run:

```bash
python tools/exclusivity_cuts.py run \
  --input /path/to/prepare_root_files \
  --diagnostics-config CLAS12/JSON_configs/rgk_epip_diagnostics.json \
  --binning-config CLAS12/JSON_configs/rgk_epip_binning.json \
  --recipe CLAS12/JSON_configs/rgk_epip_cuts.json \
  --output /path/to/exclusivity_output --workers 4
```

The stages `cache`, `derive`, `apply`, and `report` can be run separately. Model
choices, fit ranges, quality gates and window policies are configured in JSON.
Fits use only training runs. Every attempted fit and rejected/unsupported cell
is saved; each role's candidate is applied independently and exported as ROOT
friend masks. No input correction is applied. These are experimental candidate
windows, with coverage and physics-validation limits made explicit in the report.

Run all Python tests in the cuts environment with `python -m unittest discover
-s tests`. `tests/test_peakfit.py` exercises the fitter on synthetic histograms
without reading ROOT or knowing particle/sector/run identities.

## Design Philosophy

- Each particle is represented by a single configurable object (e.g. `MomCorrParticle`)
- Corrections are applied in a defined order (energy → momentum)
- Dataset- and region-dependent behavior is handled via configuration rather than hard-coding
- ROOTDataFrame supports C++ analysis; the optional Python workflow handles
  diagnostics, reusable fitting, and cut application.

## Status

Early development. APIs and structure may change significantly.  
Documentation and usage examples will be added once the core design stabilizes.
