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

## Design Philosophy

- Each particle is represented by a single configurable object (e.g. `MomCorrParticle`)
- Corrections are applied in a defined order (energy → momentum)
- Dataset- and region-dependent behavior is handled via configuration rather than hard-coding
- `ROOTDataFrame` is used throughout for scalable and transparent analysis

## Status

Early development. APIs and structure may change significantly.  
Documentation and usage examples will be added once the core design stabilizes.

