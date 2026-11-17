# TSScall benchmark harness

Submits a matrix of SLURM jobs to compare three versions of `TSScall.py` on real
data — **original** (pristine Lavender), **refac** (byte-identical bin_winner
refactor, changes 1–4), and **improved** (the full changes 1–8 build with the
terminal-uTSS fixes) — to characterize how time and peak memory scale with the
calling threshold and with input pre-filtering, and to analyze how the
called-peak universe changes across those levels.

Each version is run **bare**: the harness invokes `python <TSScall.py> …`
directly on the input bedGraphs. This does **not** go through GeneAnnotationScripts
(GGA) or any other wrapper — GGA would only matter for an integration test, which
this harness does not do. The three versions are supplied as three plain files
(see *Getting the three versions* below); you do not need three repo checkouts.

Two independent read cutoffs are swept, named consistently everywhere:

- **preFilt** — pre-filter cutoff applied to the bedGraphs *before* TSScall runs.
  A position is kept if reads ≥ preFilt. `preFilt=1` means no pre-filtering.
- **callThresh** — value passed to TSScall's `--set_read_threshold` (the internal
  calling threshold), applied to whatever input survived pre-filtering.

It answers:

1. **Do original and refac agree?** For each `(preFilt, callThresh)` cell, the
   `original` and `refac` BED body and detail file are compared by md5 — they are
   meant to be byte-identical. `improved` is **excluded from the identity check by
   design**: it differs from the other two only by the recovered terminal uTSS, so
   its delta is characterized with `analyze_overlap.py` instead.
2. **Time and memory across calling thresholds** — the `callThresh` axis, for all
   three versions.
3. **With vs. without pre-filtering** — the `preFilt` axis (`preFilt=1` is the
   "without" baseline).
4. **How the peak universe changes** — `analyze_overlap.py` reports exact /
   shifted / unique peaks across levels, split by TSS type.

## Files

| file | role |
|------|------|
| `config.sh`         | all paths, sweeps, SLURM resources, env setup. **Edit this.** |
| `submit_all.sh`     | builds pre-filtered inputs and submits the whole pipeline |
| `prefilter.sh`      | one SLURM job: strips bedGraph rows below `preFilt` |
| `run_one.sh`        | one SLURM job: runs a single condition, records time + MaxRSS |
| `collate.sh`        | one SLURM job: waits for complete metrics, then builds tables |
| `run_overlap.sh`    | one SLURM job: runs the overlap analysis after collate |
| `analyze_overlap.py`| peak-universe analytics (invoked by `run_overlap.sh`) |
| `equivalence/`      | standalone fuzz tests proving the rewritten hot paths match upstream |
| `RESULTS.md`        | distilled headline numbers from the validation run |

## Getting the three versions

The three versions are published as git **tags** on `TSScall` `main`
(`original`, `refac`, `improved`), so you extract them as three plain files into
one directory — no extra repo checkouts, no worktrees, no branch-per-version:

```bash
mkdir -p ~/tsscall_bench/versions
cd /path/to/TSScall            # your existing checkout; it stays on its branch
git fetch origin --tags
git show original:TSScall.py > ~/tsscall_bench/versions/TSScall_original.py
git show refac:TSScall.py    > ~/tsscall_bench/versions/TSScall_refac.py
git show improved:TSScall.py > ~/tsscall_bench/versions/TSScall_improved.py
```

Point the three `*_SCRIPT` paths in `config.sh` at those files. (To compare a
different pair later, extract whatever tags/commits you want the same way.)

## Quick start

1. Get the three versions (above), then edit `config.sh`:
   - `ORIG_SCRIPT`, `REFAC_SCRIPT`, `IMPROVED_SCRIPT` — the three extracted files
   - `FWD`, `REV`, `CHROMSIZES`, `GTF`, `EXPNAME`
   - `ENV_SETUP` — your cluster's `module load` / `conda activate` lines.
     TSScall is run **bare** here (not through GGA); the env just has to provide
     the interpreter and dependencies. `original` needs only Python; `refac` needs
     **numpy**; `improved` (and `analyze_overlap.py`) need **numpy + pandas**. The
     `gga` conda env already has all of these, so `conda activate gga` is the
     simplest choice — any env with numpy+pandas works equally well.
   - `PARTITION`, `JOB_MEM`, `JOB_TIME`
   - optionally the sweeps `PREFILT_LEVELS`, `CALLTHRESH_LEVELS`, `SCRIPTS`, and the
     overlap match window `OVERLAP_MAX_SHIFT`
2. Preview without submitting: `./submit_all.sh --dry`
3. Submit: `./submit_all.sh`

Everything lands under `results_<EXPNAME>/`. Jobs are chained: each pre-filter
job runs first; the TSScall jobs that use its output wait on it (`afterok`);
`collate` runs as the matrix winds down and **waits until every condition has
written its metrics row** before building tables; `run_overlap` runs after
collate succeeds.

## The matrix

```
conditions = SCRIPTS x PREFILT_LEVELS x CALLTHRESH_LEVELS
```

Defaults: `{orig, refac, improved} x {1, 4, 8} x {4, 6, 8, 10}` = 36 jobs, plus 3
pre-filter jobs, 1 collate job, and 1 overlap job. (`SCRIPTS` in `config.sh`
selects the versions; drop one to shrink the matrix.) The fixed flags mirror the
production invocation: `--annotation_join_distance 500
--annotation_search_window 1000`, with `--annotation_file`, `--detail_file`,
`--set_read_threshold`, and the input/output paths supplied per job.

Expected at full depth: the `preFilt=1` (no pre-filter) cells for **orig** OOM on
a typical node ceiling — that is the headline result, recorded as a non-zero
return code, not a harness failure — while `refac` and `improved` complete.

Note: when `callThresh < preFilt`, the pre-filter floor masks the calling
threshold (e.g. every `preFilt=8, callThresh<8` cell is identical), so only
`callThresh ≥ preFilt` cells are a clean threshold-matched comparison.

## Outputs

Under `results_<EXPNAME>/`:

- `metrics.tsv` — raw per-job row: script, preFilt, callThresh, return code,
  wall seconds, MaxRSS (`/usr/bin/time -v`), sacct MaxRSS/Elapsed cross-check,
  BED line count, BED/detail md5s.
- `comparison.tsv` — one row per `(preFilt, callThresh)`: `bed_identical`,
  `detail_identical` (YES/NO), plus orig/refac wall time and peak memory. This is
  the **orig-vs-refac identity** table; `improved` is not in it (it is expected to
  differ by the terminal uTSS — quantify that with `analyze_overlap.py`).
- `benchmark.tsv` — tidy time + peak-memory table (MaxRSS also in GB) for **all
  three** versions; the one to plot.
- `overlap_summary.tsv` / `overlap_by_type.tsv` / `overlap_shift_histogram.tsv`
  — peak-universe analytics (see below).
- `peak_dumps/<comparison>/` — `lost_from_reference.bed`, `gained_in_other.bed`,
  `shifted_reference_side.bed` for each preFilt comparison; load in IGV.
- `conditions/<script>_preFilt<PF>_callThresh<CT>/` — the BED, detail, and
  `time.txt` for each run.
- `diffs/preFilt<PF>_callThresh<CT>/` — full `bed.diff` / `detail.diff`, written
  **only** for cells that did not match (empty `diffs/` = everything agreed).
- `logs/` — per-job stdout/stderr.

## Peak-overlap analytics (`analyze_overlap.py`)

Runs automatically as the final pipeline step, and can be re-run by hand:

```
CONFIG=./config.sh bash run_overlap.sh
# or directly:
python analyze_overlap.py --results-dir results_<EXPNAME> --expname <EXPNAME> \
    --script refac --max-shift 100 --dump-dir results_<EXPNAME>/peak_dumps
```

For any pair of BED outputs it reports how many called TSSs are at identical
coordinates, how many are the same peak shifted by a few bp (with the shift
distribution), and how many are unique to one side — split by TSS type
(`obsTSS` = annotated, `uTSS` = unannotated / enhancer-class). Three comparison
families: **preFilt effect** (reference preFilt vs each other level at matched
callThresh — the read for "is no pre-filtering a reasonable standard?");
**callThresh effect** (adjacent thresholds at fixed preFilt); and **orig vs
refac** coordinate identity at every cell.

Key options: `--max-shift` sets how far apart two peaks can be and still count as
"the same, shifted" (the shift histogram always shows the true structure
regardless); `--ref-prefilt` picks the reference level (default: the lowest
present); `--dump-dir` writes the lost/gained/shifted BEDs. Needs numpy; no
bedtools required.

## Equivalence tests (`equivalence/`)

Five standalone fuzz tests prove each rewritten hot path returns exactly what the
original (or intended) code returns, on tens of thousands of randomized inputs.
They need only Python + numpy — no SLURM, no cluster, no data files — and are the
machine-checkable core of the "byte-identical" claim.

```bash
cd equivalence
python binwinner_equiv.py     # bin_winner two-pointer vs. original nested loop
python filter_equiv.py        # grouped searchsorted filter vs. intended behavior
python parse_equiv.py         # pandas/numpy bedGraph parse vs. per-line Python read
python intersect_equiv.py     # searchsorted window intersection vs. cursor walk
python int32_equiv.py         # int32 coverage arrays + an explicit >2^31 overflow check
```

Each prints a trial count and a mismatch count and ends in `PASS` (0 mismatches)
or `FAIL`. Several also run a plausible-but-wrong variant as a **negative control**,
reporting how often it diverges, to confirm the test can actually detect a bad
implementation: `binwinner_equiv` runs a fixed-tiling model, `intersect_equiv` an
overlapping-window case, and `filter_equiv` reports how often the removed filter
tail-drop would have changed the result. Expected (0 mismatches each):
`binwinner` 60,000 (tiling diverges on ~19,800), `filter` 40,000 (tail-drop fired
on ~29,700), `parse` 20,000, `intersect` 40,000 (overlapping control ~7,400),
`int32` 20,000. See [`RESULTS.md`](RESULTS.md) for the recorded numbers.

## Notes

- **The collate timing guard.** `collate.sh` counts expected conditions from the
  config and waits (up to `COLLATE_MAX_WAIT` seconds) for `metrics.tsv` to be
  complete before building tables, so a slow straggler can never yield a partial
  comparison. If it ever times out it warns loudly and you can just re-run
  `CONFIG=./config.sh bash collate.sh` once the last job finishes.
- **Memory capture.** `run_one.sh` prefers `/usr/bin/time -v` for MaxRSS; if it
  is not on `PATH`, it falls back to SLURM `sacct`. Both are recorded.
- **Pre-filtering is lossy under `bin_winner`.** Expect `preFilt > 1` cells to
  differ from the `preFilt = 1` baseline even for the same script — that is the
  effect being measured. The orig-vs-refac comparison holds preFilt fixed.
- **Threads.** TSScall is single-threaded; jobs request `--cpus-per-task=1`.
- Re-running is safe: inputs are regenerated and each condition writes to its own
  directory. Delete `results_<EXPNAME>/` for a clean slate.
