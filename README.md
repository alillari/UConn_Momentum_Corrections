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

## Design Philosophy

- Each particle is represented by a single configurable object (e.g. `MomCorrParticle`)
- Corrections are applied in a defined order (energy → momentum)
- Dataset- and region-dependent behavior is handled via configuration rather than hard-coding
- `ROOTDataFrame` is used throughout for scalable and transparent analysis

## Status

Early development. APIs and structure may change significantly.  
Documentation and usage examples will be added once the core design stabilizes.
