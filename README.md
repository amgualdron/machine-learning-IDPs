# IDP Simulation & ML

Coarse-grained molecular dynamics simulations of intrinsically disordered proteins (IDPs) using HOOMD-blue, with bulk data generation for downstream machine learning.

> **Status:** in development

---

## Overview

This project runs large batches of IDP simulations across many sequences, analyses the resulting trajectories, and builds ML-ready datasets from the outputs. Configs drive everything — any run can be fully reproduced from its YAML files saved in the run config snapshot, or the run metadata.json.

---

## Structure

```
project/
│
├── configs/                          # CONTROL PANEL
│   ├── sequences.yaml                # Full sequence database: ID -> AA string + metadata
│   ├── physics.yaml                  # Named parameter sets (default, cold)
│   └── experiment.yaml               # control panel, choose sequences to run, steps, simulation parameters
                                      # physical parameters sets or override specific parameters, choose runner
├── runs/                             # EXECUTION SPACE — fully auto-generated
│   └── exp_01_baseline/              # One folder per experiment batch
│       ├── configs_snapshot/         # Frozen copy of configs/ at generation time
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
│   └── analysis.py                   # Rg, contacts, asphericity, distributions writes run_metadata.json
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
To run a simulation locally select local as the runner in run_config.yaml, you can also specify the number
of workers (cores) that will work in parallalel for large runs, then just run the scripts normally

01*generate_jobs.py -> 02*...

it will automatically create a new directory with the run data and append entries to the master_index.csv
with all the analysis like Rg or mean Rg, the input will also be apended for ML training
**3. Submit a full experiment on the cluster**
for job summiting change the runner from local to slurm

**4. Analyse outputs**

data processing is done automatically and stored both in run_metadata.json or in the master_index.csv
, in order to generate plots or others run a custom plots.py

---

## How a run is tracked

Every simulation writes a `run_metadata.json` alongside its trajectory, for example:

```json
{
  "status": "completed",

  "exp_id": "exp_01_baseline",
  "seq_name": "ACTR",
  "param_set": "baseline",
  "run_dir": "runs/exp_01_baseline/ACTR_T300_HPS1",

  "sequence": "GTQNRPLLRNSLDDLVGPPSNLEGQSDERALLDQLHTLLSNTDATGLEEIDRALGIPELVNQGQALEPKQD",
  "length": 71,

  "temperature_K": 300,
  "hydropathy_scale": "HPS1",

  .
  .
  .

  "total_steps": 2559066,
  "n_frames": 1279,

  "wall_time_s": 3847,
  "completed_at": "2026-03-19T14:32:07"
}
```

This links every output file back to the exact code, config, and sequence that produced it.

---

## Configuration

every run is fully defined by 1 yaml file:

| File                      | Controls                                       |
| ------------------------- | ---------------------------------------------- |
| `configs/run_config.yaml` | timestep, temperature, integrator, steps, seed |
|                           | which sequences, which sim config, output dir  |

---

## Data

Raw simulation data is not tracked in git. Sequence files and all configs are. See [docs/data_management.md](docs/data_management.md) for archiving conventions.

---

## Docs

- [Data management](docs/data_management.md)
- [Cluster guide](docs/cluster_guide.md)
- [Git workflow](docs/git_workflow.md)
- [Forcefield notes](docs/forcefield_notes.md)
