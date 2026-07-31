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

## Design Philosophy

- Each particle is represented by a single configurable object (e.g. `MomCorrParticle`)
- Corrections are applied in a defined order (energy → momentum)
- Dataset- and region-dependent behavior is handled via configuration rather than hard-coding
- `ROOTDataFrame` is used throughout for scalable and transparent analysis

## Status

Early development. APIs and structure may change significantly.  
Documentation and usage examples will be added once the core design stabilizes.
