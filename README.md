# M51 VLA HI Pipeline

Scripts for reducing, imaging, and analyzing VLA HI 21-cm observations of M51.
Uses the [PHANGS imaging pipeline](https://github.com/akleroy/phangs_imaging_scripts) (`phangsPipeline`) inside CASA.

---

## Time on Source per configuration

We use the time-on-source as provided in the NRAO archive.

### 11A-142 (A-config)

These are the 11A-142 tracks we include for imaging.

| Track | Time |
|-------|------|
| 11A-142_sb4616455_1_000.55793.84399155092 | 12988.0s |
| 11A-142_sb4616455_1.55792.91121983796 | 12988.0s |
| 11A-142_sb4616455_7_002.55808.67562555555 | 12988.0s |
| 11A-142_sb4616455_7.55810.91896084491 | 12988.0s |
| 11A-142_sb4616455_1_000.55803.00196113426 | 12988.0s |
| 11A-142_sb4616455_7_000.55806.91207096065 | 12988.0s |
| 11A-142_sb4616455_7_003.55808.883936412036 | 12988.0s |
| 11A-142_sb4616455_7.55812.93744111111 | 12988.0s |
| 11A-142_sb4616455_1.55787.77720553241 | 12988.0s |
| 11A-142_sb4616455_7_001.55807.680929525464 | 12988.0s |
| 11A-142_sb4616455_7_006.55809.920632974536 | 12672.0s |

In total, there is 39.6 h of A-config time on source.


### AW605 (archival B/C/D-config)

| Track | Config | Time |
|-------|--------|------|
| AW605_1_53462.68505_53463.01734.ex | B | 6.9 h |
| AW605_1_52942.11259_52942.42174.exp | B | 6.2 h |
| AW605_1_53057.85773_53058.41676.exp | C | 2.0 h |
| AW605_1_53191.67833_53192.09245.exp | D | 1.3 h |

## Summary by Config

| Config | Time (h) |
|--------|:--------:|
| A      | 39.6     |
| B      |  13.1     |
| C      |  2.0     |
| D      |  1.3     |
| All    | 56     |

---

## Calibration

Two separate calibration workflows handle the modern (11A-142) and archival (AW605) data.

### 11A-142

**Scripts:** `11A-142/`

- `11A-142_manual_lband_script.py` — Full manual L-band calibration pipeline for a single track. Steps: Hanning smoothing, edge/quack flagging, per-track manual flags (from `manual_flags/`), initial phase calibration, delay calibration, bandpass calibration, amplitude gain calibration, fluxscale, applycal, rflag on continuum, and final target split.
- `run_pipeline_11A-142.sh` — Wrapper to process a single track (extracts tar, runs CASA, runs QA plots via `qaplotter`).
- `runner_11A-142.sh` — Runs all tracks in parallel via GNU `parallel` (4 at a time).

```bash
# Run all A-config tracks in parallel:
cat ~/M51_HI/11A-142/track_names.txt | parallel -n 4 ~/M51_HI/11A-142/run_pipeline_11A-142.sh
```

### AW605

**Scripts:** `AW605/`

- `AW605_1_*_reduction.py` — One script per track/config. Imports data from AIPS `.exp` files via `importvla`, calibrates, and splits the science target with a test dirty cube.
- `add_scan_intents.py` — Utility to patch missing scan intents into old MS files before calibration.

---

## Imaging

### HI Line Imaging

**Script:** `run_casa_pipeline_hi.py` | **Job array:** `jobarray_imaging_hi.job`

Must be run inside CASA. Controlled by three flags: `do_uvprocessing`, `do_imaging`, `do_postproc`.

**UV Processing** (`do_uvprocessing`): Copies MS files, runs continuum subtraction (`uvcontsub`), extracts line data per config, and applies `statwt` noise weighting.

**Imaging** (`do_imaging`): To avoid CASA memory limits with `tclean`, the full M51 field is split into 9 channel chunks (`m51_1`–`m51_9`). Each chunk is imaged with `loop_imaging` using recipe `phangsalma`. OH lines receive dirty images only; HI lines use full multiscale clean. Mosaicking along the spectral dimension must be done before post-processing.

Line products (SGE task IDs 1–5):

| ID | Product | Configs |
|----|---------|---------|
| 1 | `hi21cm_2p5kms` | A |
| 2 | `hi21cm_5kms` | A, B+C+D, C+D, D |
| 3 | `hi21cm_10kms` | A, B+C+D, C+D, D |
| 4 | `oh1665` | A |
| 5 | `oh1667` | A |

```bash
# Submit imaging job array (tasks 1–5):
qsub jobarray_imaging_hi.job
```

### Continuum Imaging

**Scripts:** `continuum_imaging/`

Three-stage `tclean`: single-scale clean → multiscale MS-clean of the central region → PB-masked single-scale clean, followed by `widebandpbcor`. Uses `wproject` gridder (256 planes), `mtmfs` deconvolver (`nterms=2`), Briggs robust=0.5, 0.2″/cell, 16384-pixel image.


---

## Derived Products

**Script:** `run_derived_pipeline_hi.py` (run outside CASA)

Runs the full derived-product pipeline on the mosaicked HI cubes for all interferometric configs (`D`, `C+D`, `B+C+D`, `A+B+C+D`) and HI line products (`hi21cm_5kms`, `hi21cm_10kms`).

---

## Analysis Scripts

These scripts run in a standard Python environment (`spectral-cube`, `astropy`, `regions`).

| Script | Description |
|--------|-------------|
| `measure_HI_flux.py` | Measures integrated HI flux density and HI mass from moment-0 maps; checks flux recovery vs. archival literature measurements |
| `measure_HI_spectrum.py` | Extracts and plots spatially integrated HI spectra for flux/mass verification |
| `plot_hi_summary.py` | Produces summary diagnostic plots of the HI cubes and moment maps |
