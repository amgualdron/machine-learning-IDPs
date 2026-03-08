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
├── data/
│   ├── sequences/          # Sequence files (.dat) + metadata.yaml
│   ├── raw/                # HOOMD outputs — trajectory.gsd, thermo.log, run_manifest.yaml
│   └── processed-ml/       # Feature matrices ready for ML training
│
├── configs/
│   ├── simulation/         # Integrator, timestep, temperature, steps
│   └── experiments/        # Which sequences to run + which sim config to use
│
├── src/
│   ├── config.py           # YAML loading, sequence translation, validation
│   ├── simulation.py       # HOOMD runner
│   ├── analysis.py         # Trajectory observables (Rg, contact maps, ...)
│   ├── features.py         # Processed data → ML feature matrices
│   └── train.py            # ML model training
│
├── scripts/
│   ├── run_simulation.py   # --experiment configs/... --index N
│   ├── run_analysis.py
│   └── build_features.py
│
├── slurm/
│   ├── submit_experiment.sh   # SLURM job array → calls run_simulation.py
│   └── submit_analysis.sh
│
├── logs/                   # Cluster stdout/stderr, gitignored
├── notebooks/              # Exploratory analysis only, not part of pipeline
├── tests/
└── docs/
└── enviroment.yaml # enviroment packages to recreated using conda
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
