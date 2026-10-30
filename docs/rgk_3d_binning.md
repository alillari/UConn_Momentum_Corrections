# RG-K sector-separated 3D binning

The first step toward exclusivity fits is a cell scheme in **momentum magnitude,
theta, and shifted phi**, defined independently for electrons and positive
pions. Sector labels are separate categories and are never pooled in the saved
calibration cells. No exclusivity cut, fit, or correction is applied here.

## Coordinates and intervals

Momentum is in GeV/c and angles are in degrees. Azimuth comes from `atan2(py,px)`.
Following the existing convention, add 360 degrees when sectors 3/4 have a
negative azimuth, or sectors 5/6 have azimuth below 90 degrees. Subtract
`(sector - 1)*60` to obtain the sector-local azimuth.

- Electron: `phi_shift = phi_local - 30/p`.
- Positive pion: `phi_shift = phi_local + 32/(p - 0.05)`.

These equations were confirmed by the user. They correspond to the electron
and forward-pion expressions in [Richard's source](https://github.com/RichCap/Exclusive_RG-A_Momentum_Corrections/blob/5815fa289f6a9bbae3c2887ccabb364edb3f2e85/Data_Files/Momentum_Correction_File_Creation_wPass2.py#L1505).
Their use here follows the user's chosen coordinate convention for the supplied
sector labels; a detector-region classification is not inferred from this skim.
There is no wrapping after the inverse-momentum shift. Invalid sector labels,
zero/nonfinite momenta, and momenta at/below the shift pole are explicit invalid
categories. Bins are half-open with the final right edge included.

## Initial selected grids

The intermediate (`balanced`) schemes retain useful resolution in populated
regions and widen theta/momentum bins in sparse tails. They are starting
calibration grids, not proof that every cell supports a reliable fit.

| Coordinate | Electron | Positive pion |
| --- | --- | --- |
| Momentum | 0.5-6.75, width 0.25 | 0.25-3, width 0.25; then edges 3.5, 4, 5, 6.5 |
| Theta edges | 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 24, 28, 32, 40 | 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 52, 60, 75, 90 |
| Shifted phi | -30 to 30, width 10 | -40 to 40, width 10 |
| Cell shape | 6 sectors x 25 p x 13 theta x 6 phi | 6 sectors x 15 p x 14 theta x 8 phi |

The exact edge vectors and the coarse/fine alternatives are versioned in
`CLAS12/JSON_configs/rgk_epip_binning.json` and copied into each output bundle.

## Full-data comparison

All 471 files and 15,731,576 entries were used. A **reference mass band of
0.85-1.05 GeV** is counted only to assess local population near the neutron
region. It contains background and is not applied as a selection. No
signal-yield, purity, efficiency, or fit-quality claim follows from it.

| Particle / scheme | Occupied cells | Reference-band entries in cells with >=300 band entries |
| --- | ---: | ---: |
| Electron / coarse | 1,941 | 96.76% |
| Electron / intermediate | 7,531 | 90.51% |
| Electron / fine | 27,713 | 73.70% |
| Pion / coarse | 1,237 | 98.54% |
| Pion / intermediate | 4,869 | 93.92% |
| Pion / fine | 28,485 | 35.46% |

Percentages are conditional on lying inside the respective grid. The selected
electron grid covers 99.99922% of input entries; the pion grid covers 99.92304%.
At a 100-entry reference-band threshold, supported cells contain 95.78% and
98.26% of their in-grid reference entries, respectively. The fine pion grid has
no cell with 1,000 reference-band entries, so uniformly fine binning would make
the initial fits unnecessarily sparse. Coarse alternatives remain available
for explicit, reviewable changes in sparse regions; there is no silent merge
or extrapolation across cells.

Outside the selected electron grid are 4 momentum-domain, 70 theta-domain, and
48 shifted-phi-domain entries. Outside the pion grid are 140, 81, and 11,886,
respectively. These are disjoint categories, assigned in the order p, theta,
shifted phi. They are not discarded from the input accounting.

The input contains rare extreme pion momenta, reaching about 399 GeV/c. This
is outside the beam-energy scale and calibration domain. The summary records
file/run/event examples at coordinate extrema for tracing these rows. Their
cause is not established by the ten available branches.

## Export and verification

Each selected-grid cell has a missing-mass spectrum from 0.65 to 1.35 GeV in
2 MeV bins, plus an explicit invalid/out-of-range count. Spectrum counts plus
those losses equal its occupancy. Counts by run are retained for subsequent
run-separated validation. The saved CSVs include empty cells and all exact
boundaries; the NPZ stores all three candidate population grids and the
selected-grid spectra. Report projections sum sectors for display only.

The tests cover shifted-phi signs/sector unwrap, degree units, singular shifts,
all three dimensions' edge rules, separate sector IDs, out-of-domain accounting,
strict schema validation, and ROOT-to-cell/spectrum integration. All eleven
Python tests pass, including the earlier diagnostic physics checks.

Next, assess peak/background fits on these individual cell spectra. Sparse or
ambiguous cells must be reported explicitly rather than receiving a cut from
an unrelated projection or neighboring sector.
