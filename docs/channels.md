# Developing cuts for charged particles and channels

The schema-2 cut workflow distinguishes **electron, pi_plus, pi_minus, and
proton**. Particle identity, measured topology, calibration role, coordinate
grid, and fitted observable are separate declarations. The same cache,
statistical fitter, window table, mask exporter, and report serve each channel.

| Example directory | Final state | Measured tracks | Example observable |
| --- | --- | --- | --- |
| `examples/channels/epipn` | e pi+ n | e, pi+ | Missing mass of n |
| `examples/channels/ep_elastic` | e p | e, p | Full missing energy |
| `examples/channels/eppippim` | e p pi+ pi- | e, p, pi+, pi- | Full missing energy |

**pi_minus is accepted only in the fully measured ep pi+ pi- channel.** An
absent pi- branch is an error; it cannot become an inferred measured track.
A nonfinite or zero declared measured momentum invalidates the observable,
including in leave-one-out calculations. Each input row must already be a
selected candidate with the declared particle identities. This adapter does
not perform PID, check candidate multiplicity, or reconstruct HIPO tracks.
Branch mappings cannot alias pi+ and pi- to the same source fields.

These examples are **templates, not validated cuts**. Only epipn is available
locally as real data. Elastic and four-track channels have been checked with
conservation-law fixtures and synthetic peak/background samples, not measured
calibration data. Match their example branch names to the actual skims.
Conditions, PID/QA history, corrections, run partitions, spectra and bin edges
require confirmation for each dataset.

## What a new cut changes

1. `input.json` declares final state, measured tracks, logical branch names,
   beam total energy, proton-at-rest target, and explicit particle masses.
2. `binning.json` declares each calibration particle's sector-separated
   p/theta/shifted-phi grid and its own shift coefficient and pole.
3. `cuts.json` chooses the observable, spectrum/units, calibration roles,
   training/validation runs, fit family, quality gates, and window rule.

For example, to study the fully measured elastic channel:

```bash
python tools/exclusivity_cuts.py run \
  --input /path/to/elastic_root_files \
  --dataset-config examples/channels/ep_elastic/input.json \
  --binning-config examples/channels/ep_elastic/binning.json \
  --recipe examples/channels/ep_elastic/cuts.json \
  --output /path/to/elastic_candidates \
  --workers 4
```

For the fully measured double-pion channel, use the three files under
`examples/channels/eppippim/` and a corresponding four-track skim. The current
two-track RG-K files cannot supply that channel or a measured elastic proton.
`--dataset-config` and `--diagnostics-config` are aliases; the former is clearer
for schema 2. Original schema-1 commands and saved artifacts remain supported,
with `pion` explicitly meaning the legacy positive pion. Schema 2 requires
charge-specific identities. Mixing schemas across input, grids and recipe is
rejected.

## Observables and topology

For an electron beam along +z and a proton target at rest, define

```
R = beam + target - sum(included measured four-vectors).
```

Every track uses its declared species mass for energy. Available kinds are:

| Kind | Value | Unit |
| --- | --- | --- |
| `missing_mass` | signed sqrt(R.E^2 - R.p^2) | GeV/c^2 |
| `missing_mass2` | R.E^2 - R.p^2 | GeV^2/c^4 |
| `missing_energy` | R.E | GeV |
| `missing_px`, `missing_py`, `missing_pz` | signed residual component | GeV/c |
| `missing_p` | magnitude of R.p | GeV/c |
| `missing_pt` | transverse magnitude of R.p | GeV/c |

With all final tracks measured, ideal full missing energy and momentum close
at zero. They do not have the neutron reference used for epipn. The fit family
must suit the observable: a nonnegative magnitude or signed square-root mass
near zero is not automatically Gaussian when detector components are Gaussian.

Choose full energy closure with:

```json
{"kind": "missing_energy", "omit": []}
```

`omit` excludes measured tracks **from the four-vector subtraction**. For a
leave-one-out proton mass in elastic data:

```json
{"kind": "missing_mass", "omit": ["proton"]}
```

For a leave-one-out negative-pion mass in the measured four-track channel:

```json
{"kind": "missing_mass", "omit": ["pi_minus"]}
```

The omitted track is still required in the input. The latter is
M(beam + target - e - p - pi+), studied versus measured pi- coordinates;
it does not enable a missing-pi- input topology. A leave-one-out mass does not
directly test the omitted track's measured four-vector. Correlation with its
calibration coordinates and other selections must be considered. An omitted
electron's tiny mass makes a signed square-root extraction numerically and
experimentally delicate; mass squared preserves the algebraic constraint.

`spectrum` declares edges, a label, units, and an optional reference line.
Reference values are annotations, not fit constraints. The center is free
within its configured range. Changing the observable, omissions, spectrum,
run split, or roles requires a new cache. Fit/background/window changes can
reuse it. One workflow fits one observable; run separate recipes for others.

## Coordinates and present boundaries

All current grids are explicitly **FD**, with sector IDs 1 through 6. CD is
rejected rather than passed through the FD azimuth convention. CD coordinates
and mixed FD/CD calibrations need a separate grid adapter.

Electron and pi+ examples retain the approved shifts: local phi - 30/p and
local phi + 32/(p - 0.05), in degrees and GeV/c. There is no charge-based guess
for proton or pi- shifts. Their templates explicitly use zero coefficient
(unshifted sector-local phi) as an untuned starting point. Appropriate grids
and conventions must be established from their own data. Coverage, fit support,
and cut rejection are reported separately.

Masks remain independent per calibration role. These are selections conditional
on each particle's coordinates, not independent exclusivity measurements.
Applying an explicit AND/OR across roles or observable selections is a separate
analysis decision. Automatic cut composition, joint efficiency optimization,
and ML correction training are not implemented.

## Implementation and verification

`momcorr.channels` owns topology checks and four-vector observables. Binning,
caching, evaluation and reporting read the declared roles and spectrum. The
fitter and window table remain independent of physics channels. Reports split
particle comparisons into groups of at most two.

Tests generate exact elastic and sequential double-pion decays and check
four-vector closure, omitted-track masses, charge/identity constraints,
missing/invalid tracks, independent shifts, cache compatibility, and aligned
per-particle masks. Synthetic residual peaks plus background pass through
caching, fitting and application for elastic and four-track inputs. This tests
the software contract and recovery, not physics validity on measured data.

The complete suite has 34 passing tests. Independent CLI runs with two fitting
workers produced masks and reports for 9,000 synthetic rows per new channel.
Rebuilding the canonical epipn cache over all 471 real files (15,731,576 rows)
exactly reproduced the original electron/pi+ training and validation spectra,
grid populations, and out-of-spectrum populations. Existing schema-1 results
also render with the updated report code. No new real-data fit or validated
elastic/double-pion correction is implied by these checks.
