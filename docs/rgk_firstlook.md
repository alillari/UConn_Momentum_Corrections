# RG-K single-pion first diagnostic pass

This is the first input and kinematic audit for the exclusivity-cut study. It
does not establish an exclusivity selection or derive a momentum correction.

## Input and conditions

The local input comprises 471 ROOT files, 410,196,038 bytes, and 15,731,576
entries in one consistent `h22` schema. All `(run, evn)` pairs are unique.

| Run | Files | Entries |
| --- | ---: | ---: |
| 5875 | 33 | 1,355,444 |
| 5877 | 13 | 627,523 |
| 5880 | 36 | 1,392,795 |
| 5896 | 197 | 6,270,557 |
| 5932 | 192 | 6,085,257 |

The user identified these as RG-K fall 2018, 6535 MeV, Pass 2, and confirmed
that the supplied momenta have no energy-loss corrections. Other preparation
selections and correction history are not encoded in these ROOT files.

The nominal total beam energy is **6.53536 GeV** and the target is LH2, modeled
as a proton at rest. The source is the [JeffersonLab RG-K configuration](https://github.com/JeffersonLab/clas12-config#run-configurations).
The [RG-K analysis note](https://clasweb.jlab.org/wiki/images/e/ed/Rgk_analysis_document.pdf)
maps runs 5875-6000 to the 6.535 GeV period and documents outbending production
conditions. This note describes pass-1 processing: its data-taking conditions
are relevant, but its cuts and QA certification are not adopted for this
pass-2 sample. Per-run RCDB beam energies were not retrieved.

## Observations

- A peak consistent with the neutron and its band in missing mass squared are
  visible in the inclusive sample and all six electron/pion sector groups.
- The higher-mass continuum is substantial. A peak window needs a local
  background study and should not be chosen from a single global Gaussian fit.
- The electron low-momentum turn-on differs across runs. Run-separated
  comparisons must be retained while the preparation history is investigated.
- No nonfinite momenta or zero-momentum candidates were found. There are
  15,129 negative missing-mass-squared entries, 3,519 negative missing energies,
  and 27 electron momenta above the nominal beam energy. These rows remain in
  the analysis; display overflows are recorded rather than silently discarded.
- All input sector labels are in 1-6. Detector status, PID, vertex, timing,
  track quality, and correction-history branches are absent. Sector labels
  alone are not used to classify detector regions.

The pion on-shell-neutron constraint has 15,714,146 unique physical roots,
9,595 ambiguous cases, and 7,835 cases with no accepted physical root. This is
an estimator under a channel hypothesis, not truth. Only the unique-root cases
enter its delta-p plot; no momentum is replaced and no other plot is selected
on this result. Candidate roots are checked against the unsquared energy
equation, not selected for proximity to the measured pion momentum.

## Validation

- Four automated tests pass, covering a generated two-body final state,
  preservation of negative mass squared, complete histogram accounting, and
  ROOT integration with repeated candidate IDs.
- On 1,000 actual entries, the NumPy signed missing mass agrees with the
  existing `include/Kinematics.h` calculation compiled with C++17 to a maximum
  absolute difference of approximately `5.33e-15 GeV`.
- All 65 histogram input/display/overflow/nonfinite identities balance.
  Manifest and run totals agree with the entry total.
- The diagnostic PDF is visually reviewed. The original input files are only
  read, and their SHA-256 hashes are saved in the output manifest.

## Next step

Fit the neutron band in local kinematic/sector bins, first with a small set of
reviewable peak/background alternatives. Compare candidate boundaries on
held-out runs and measure their effect on the pion momentum-residual spectrum.
No cut values or efficiency/purity claims are established in this first pass.

The reproducible command and dependency list are in the repository README.
Report, histograms, input hashes, code/configuration snapshot, and display
accounting are generated outside the source tree.
