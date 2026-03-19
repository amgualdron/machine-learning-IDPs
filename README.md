# IDP Simulation & ML

Coarse-grained molecular dynamics simulations of intrinsically disordered proteins (IDPs) using HOOMD-blue, with bulk data generation for downstream machine learning.

> **Status:** in development — physics lab collaboration

---

## Overview

This project runs large batches of IDP simulations across many sequences, analyses the resulting trajectories, and builds ML-ready datasets from the outputs. Configs drive everything — any run can be fully reproduced from its YAML files and a git commit hash.

---

## Structure

```
project/
│
├── configs/                          # MASTER CONTROL PANEL — hand-edit here only
│   ├── sequences.yaml                # Full sequence database: ID -> AA string + metadata
│   ├── physics.yaml                  # Named parameter sets (baseline, cold, ph_sweep...)
│   └── experiment.yaml               # Which sequences x param_sets to run + sweep flags
│
├── runs/                             # EXECUTION SPACE — fully auto-generated
│   └── exp_01_baseline/              # One folder per experiment batch
│       ├── configs_snapshot/         # Frozen copy of configs/ at generation time
│       │   ├── sequences.yaml        #   (full DB, not just sequences that ran)
│       │   ├── physics.yaml
│       │   └── experiment.yaml
│       ├── sequences_ran.yaml        # Explicit list of sequences in THIS experiment
│       ├── manifest.csv              # One row per job — used by SLURM array indexing
│       └── fus_lcd_T300_HPS1/        # Auto-named: {seq}_{T}K_{hps_scale}
│           ├── run_metadata.json     # Written TWICE: pre-run (config) + post-run (status)
│           ├── trajectory.gsd        # Raw HOOMD output
│           └── hoomd.log             # Raw HOOMD output
│
├── data/
│   ├── external/                     # Downloaded files, PDBs, reference data
│   └── processed-ml/                 # Distilled, ML-ready outputs
│       ├── exp_01_features.csv       # Computed observables (Rg, contacts, etc.)
│       └── exp_01_targets.npy        # Numpy arrays for large per-frame data
│
├── scripts/                          # Workflow automation — numbered in order
│   ├── 01_generate_jobs.py           # configs/ → runs/exp_N/ + manifest.csv
│   └── 02_process_to_ml.py           # runs/ → data/processed-ml/ + master_index.csv
│
├── src/                              # Core library — importable by scripts
│   ├── simulation.py                 # HOOMD script, reads run_metadata.json
│   ├── analysis.py                   # Rg, contacts, asphericity, distributions
│   └── config_loader.py              # YAML loading, deep merge, validation, extends
│
├── slurm/
│   └── submit_array.sh               # Indexes into manifest.csv via $SLURM_ARRAY_TASK_ID
│
├── master_index.csv                  # ONE ROW PER COMPLETED+ANALYZED RUN (see below)
├── environment.yml
└── README.md
```

---

## Quickstart

**1. Set up the environment**

```bash
git clone https://github.com/<you>/project.git
cd project
conda env create -f environment.yml
conda activate idp-research
```

**2. Run a simulation locally**

```bash
python scripts/run_simulation.py \
    --experiment configs/experiments/charge_sweep.yaml \
    --index 0
```

**3. Submit a full experiment on the cluster**

```bash
sbatch slurm/submit_experiment.sh configs/experiments/charge_sweep.yaml
```

**4. Analyse outputs**

```bash
python scripts/run_analysis.py --experiment configs/experiments/charge_sweep.yaml
```

**5. Build ML features**

```bash
python scripts/build_features.py \
    --processed-dir data/processed-ml/ \
    --out-dir data/processed-ml/features/
```

---

## How a run is tracked

Every simulation writes a `run_manifest.yaml` alongside its trajectory:

```yaml
run_id: charge_sweep_polyR10_20260306_143201
git_commit: a3f9c2d
experiment_config: configs/experiments/charge_sweep.yaml
sequence_file: data/sequences/polyR10.dat
sequence_hash: md5:e4d909c290d0fb1ca068ffad
slurm_job_id: 10042
timestamp: 2026-03-06T14:32:01
```

This links every output file back to the exact code, config, and sequence that produced it.

---

## Configuration

Two YAML files fully define any run:

| File                              | Controls                                       |
| --------------------------------- | ---------------------------------------------- |
| `configs/simulation/default.yaml` | timestep, temperature, integrator, steps, seed |
| `configs/experiments/<name>.yaml` | which sequences, which sim config, output dir  |

---

## Data

Raw simulation data is not tracked in git. Sequence files and all configs are. See [docs/data_management.md](docs/data_management.md) for archiving conventions.

---

## Docs

- [Data management](docs/data_management.md)
- [Cluster guide](docs/cluster_guide.md)
- [Git workflow](docs/git_workflow.md)
- [Forcefield notes](docs/forcefield_notes.md)
