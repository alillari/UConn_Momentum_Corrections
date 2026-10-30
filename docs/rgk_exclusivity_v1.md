# RG-K first candidate exclusivity windows

This exploratory study used the current **471-file, 15,731,576-entry snapshot** and source commit `ee8dde8c1d98d45e4bc1098f4f25490fb63fb8c3`. It establishes a working, reusable fit-to-window-to-mask pipeline; these candidates are not production-validated physics cuts.

## Data and evaluation scope

The supplied file list and local directory both contain 471 ROOT files. The reported expectation of 500 files has not been reconciled upstream. The user authorized continuing with the available sample. All input hashes remain unchanged.

Training runs: 5875, 5877, 5896. Runs withheld from fitting: 5880, 5932. Binning was explored on the full sample before this split; this is exploratory run-held-out evaluation, not a blind final test.

The supplied preparation script indicates electron identification and exactly one accepted pion, with QA disabled in that script. External PID definitions and identity with the actual producer code are unverified. The absence of energy-loss corrections is user-confirmed. Input momenta are preserved.

## Candidate results

| Selection role | Ready cells | Review cells | Insufficient cells | Validation rows in supported cells | All validation rows retained |
| --- | ---: | ---: | ---: | ---: | ---: |
| electron | 1,245 | 130 | 10,325 | 48.95% | 9.48% |
| pion | 1,734 | 184 | 8,162 | 92.97% | 13.09% |

Insufficient counts include empty cells and cells with too few training entries in the diagnostic mass band. `ready` means the configured algorithmic gates passed, not that physics validity was established. Electron and pion decisions are independent; they are not silently combined. Pass fractions are retention, not signal efficiency or purity.

The electron-based candidate has a substantial coverage limitation. Unsupported regions receive no fallback cut, interpolation, or neighboring-sector window. Expanding coverage requires an explicit new binning/data/model study.

## Method and saved evidence

The recipe uses bin-integrated Poisson Gaussian-plus-positive-background fits in three mass ranges, with polynomial degrees 0, 1 and 2. AICc comparisons stay within the same range; plausible models and range winners must have stable centers and widths. The nominal window is fitted center ±3 fitted widths. No per-cell special cases or manual fit adjustments were used.

All fit attempts and failures are saved. Review gates cover optimizer/covariance problems, weak peaks, structured residuals, boundary fits, model/range sensitivity, and cuts extending beyond the fitted domain. ROOT friend masks preserve original row order and run/event IDs; separate reason codes distinguish unsupported cells from rejected events.

Before/after missing mass is a selection diagnostic on the selected variable. The inferred pion momentum residual is correlated with missing mass, so narrowing it is not independent correction validation. Independent signal-efficiency/purity and unbiased correction claims remain open.

## Verification

- All 24 Python unit/integration tests pass, including injected peak recovery, model derivatives, unit scaling, validation isolation, and friend-mask alignment.
- Cached training-plus-validation populations and mass spectra exactly reproduce the earlier full-data 3D-binning archive.
- Full masks were checked against source entry order and run/event IDs, source hashes and per-run evaluation counts.
- In 80 correctly specified Gaussian/sloped-background toys, all fits passed and nominal 68% center intervals covered 67.5%. This does not establish coverage under arbitrary backgrounds.

The architecture, CLI stages, output formats and reuse boundaries are described in `docs/exclusivity_workflow.md`. The PDF report includes accepted/review examples, residuals, parameter distributions, held-out comparisons, and width sensitivity.
