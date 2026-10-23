# Reproduction log: final verification

## Environment

| Item | Value |
|---|---|
| Machine | Cloud container (Linux) |
| R | 4.3.3 (2024-02-29) |
| Key packages | survival 3.5.8, rms 6.7.1, glmnet 4.1.8, mice 3.16.0, data.table |
| Full session info | `audit/output/sessionInfo.txt`, `post_audit/output/sessionInfo.txt` (unchanged, same session) |
| Repository state at start | commit `21e94f0`; working tree clean; branch `claude/uterine-pathology-stroke-stats-fb93xb` |
| Interim commit during verification | `ebe420e` (final_verification files only) |
| Data | `../data` (unchanged since 2026-09-29) and `derived/*.rds` from `run_all.R` |

## Seeds and resampling

All scripts use `SEED = 20261022` (from `R/00_setup.R`).

| Script | Resampling | Seed |
|---|---|---|
| V01 | 1,000 paired resamples of the 2023+ set | `SEED` |
| V01 | 1,000 paired resamples of the baseline cohort | `SEED + 1` |
| V04 | 1,000 resamples within each fibroid subgroup | `SEED` |
| V02 | 200 patient-level resamples × 7 models (P05 logic) | `SEED` |

V02 also re-runs P05's apparent fits, which include `cv.glmnet`; `set.seed(SEED)` precedes them, as fixed in P05 on 2026-10-23.

## Commands run (from `reanalysis/`)

```
Rscript final_verification/R/V01_reconcile.R
Rscript final_verification/R/V02_rerun_p05_instrumented.R      # ~25 min; writes derived/final_verification/P05_bootstrap_rerun.rds
Rscript final_verification/R/V03_drift.R
Rscript final_verification/R/V04_fibroid.R
Rscript final_verification/R/V05_sparse_calculator.R
Rscript final_verification/R/V06_compare_bootstrap.R           # requires V02
```

All six completed with exit code 0. The last runs of V03, V04 and V06 came after the corrections listed below.

## Corrections made to my own verification code during the run (no original file affected)

- **V00:** added `require_cols()`, which was missing on first run (V01 stopped with "could not find function").
- **V03, row removed:** the "all calendar-2023+ follow-up" row for pre-2023 anchors relied on the baseline hazard beyond the development model's 5-year range (held flat), which understates expected counts. It was dropped and the reason is documented in V03.
- **V03, window fixed:** the old-cohort window was first capped at 3 years after entry, so up to 6 years after the anchor. It was corrected to 0–3 years since the anchor for both cohorts. This changed that row's O/E from 0.97 to **0.59**; only the corrected value is reported.
- **V06, note corrected:** the note on zero-event resamples first pooled the two outcomes. It was corrected to count each outcome separately.

## Warnings

- **V02 (P05 re-run):** 4 coxph warnings, "Loglik converged before variable …; coefficient may be infinite". These are the same 4 as the original P05 run. They are diagnosed in V06b as quasi-separation (0 events among `hbp_pre`-positive women in the resample).
- **All scripts:** "strings not representable in native encoding" when reading `anchored_cohort.rds`. This is cosmetic (en-dashes in labels).

## Failed or blocked steps

- **Loading `outputs/_state.pkl`:** blocked by the session's permission policy (an unpickling risk) in the post-audit stage. Its source (`analysis/run_all.py`) shows it holds result tables only, so it was not needed.
- **Live web page in the headless browser:** not reachable because of the proxy certificate (audit stage). The page was fetched with curl and compared byte for byte instead.
- **Dating of migraine, diabetes and estrogen; outcome PPV:** not attempted. The data are absent (`MISSING_DATA_REQUEST.md`).

## Preservation check

- `final_verification/checksums_before.txt` and `checksums_after.txt` hold MD5 sums of all 243 files under `output/`, `audit/output/`, `post_audit/output/`, `calculator/` and `derived/` (excluding `derived/final_verification/`). They were taken before and after all verification runs, and **are identical**.
- V02 executes a *copy* of `post_audit/R/P05_candidate_models.R` with three documented string substitutions: the cache path, the output directory, and one recorded variable. The original script is untouched. Sourcing `P00_setup.R` rewrites `post_audit/output/sessionInfo.txt` with identical content, and its checksum is unchanged.
- The re-run P05 report (`final_verification/output/p05_rerun/P05_candidate_models.md`) is byte-identical to `post_audit/output/P05_candidate_models.md`, apart from the generated-on header.

## Not reproduced in this verification (relied on as previously run)

- **Audit scripts A01–A11:** last clean run on 2026-10-23, commit `cd77c82`.
- **Post-audit P01–P04 and P06–P08:** last clean run at commit `21e94f0`.

Their key numbers were independently recomputed in V01 (temporal C, IPA, O/E, follow-up, event counts) and matched.
