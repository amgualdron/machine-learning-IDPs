import argparse
import csv
import itertools
import json
import shutil
from datetime import datetime
from pathlib import Path

import yaml


# ── Helpers ───────────────────────────────────────────────────────────────────
def load_yaml(file_path: Path) -> dict:
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


# ── Main Generator Logic ──────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate job directories and manifests."
    )
    parser.add_argument(
        "--config", required=True, help="Path to the master run_config.yaml"
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    project_root = config_path.parent.parent

    run_config = load_yaml(config_path)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    exp_name = run_config.get("run_name", "unnamed_run")
    exp_id = f"{timestamp}_{exp_name}"
    # 1. Load Sequences
    seq_source_path = config_path.parent / run_config["sequences"]["source"]
    all_sequences = load_yaml(
        seq_source_path
    )  # Assumes format: { "ACTR": "ASDF...", ... }

    # 2. Sequence Selection Logic
    selection = run_config["sequences"].get("select", "all")
    selected_seqs = {}

    if isinstance(selection, list):
        # Explicit list provided
        for name in selection:
            if name in all_sequences:
                selected_seqs[name] = all_sequences[name]
            else:
                print(f"Warning: Sequence '{name}' not found in {seq_source_path.name}")

    elif selection == "all":
        # Full database
        selected_seqs = all_sequences

    # 3. Setup Experiment Directory
    exp_dir = project_root / "runs" / exp_id
    exp_dir.mkdir(parents=True, exist_ok=True)
    snapshot_dir = exp_dir / "configs_snapshot"
    snapshot_dir.mkdir()
    shutil.copy(config_path.parent / "physics.yaml", snapshot_dir)
    shutil.copy(config_path, snapshot_dir)
    with open(snapshot_dir / "ran_sequences.yaml", "w") as f:
        yaml.dump(selected_seqs, f)

    # Create logs directory for SLURM output
    (exp_dir / "logs").mkdir(exist_ok=True)

    manifest_path = exp_dir / "manifest.csv"
    manifest_rows = []

    # 4. The Cartesian Product Loop (Sequences x Parameters)
    task_id = 1
    for seq_name, param_set in itertools.product(
        selected_seqs.keys(), run_config["param_sets"]
    ):
        seq_data = selected_seqs[seq_name]
        seq_string = seq_data["sequence"]

        # Create specific run directory
        run_name = f"{seq_name}_{param_set}"
        run_dir = exp_dir / run_name
        run_dir.mkdir(exist_ok=True)

        # Write Stage 1 run_metadata.json
        meta = {
            "experiment_id": exp_id,
            "sequence_name": seq_name,
            "sequence": seq_string,
            "param_set": param_set,  # Singular string now, so resolve() works in simulation.py
            "status": "pending",  # Ready to be updated by simulation.py
        }

        with open(run_dir / "run_metadata.json", "w") as f:
            json.dump(meta, f, indent=2)

        # Append to manifest tracking
        manifest_rows.append(
            {
                "task_id": task_id,
                "run_dir": str(run_dir),
                "sequence_name": seq_name,
                "length": len(seq_string),
                "param_set": param_set,
            }
        )
        task_id += 1

    # 5. Write the Master Manifest
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["task_id", "run_dir", "sequence_name", "length", "param_set"]
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    # 6. Generate Execution Script
    total_jobs = len(manifest_rows)
    runner = run_config.get("runner", "local")

    if runner == "slurm":
        slurm = run_config["slurm"]
        max_concurrent = slurm.get("max_concurrent", 32)

        submit_script = f"""#!/bin/bash
#SBATCH --job-name={exp_id}
#SBATCH --output={exp_dir}/logs/%A_%a.out
#SBATCH --error={exp_dir}/logs/%A_%a.err
#SBATCH --array=1-{total_jobs}%{max_concurrent}
#SBATCH --partition={slurm['partition']}
#SBATCH --gres={slurm['gres']}
#SBATCH --time={slurm['time']}
#SBATCH --mem={slurm['mem']}

# Extract run_dir from the Nth line of the manifest (adding 1 to skip the header)
ROW_NUM=$(($SLURM_ARRAY_TASK_ID + 1))
RUN_DIR=$(sed -n "${{ROW_NUM}}p" {manifest_path} | cut -d',' -f2)

echo "Starting task $SLURM_ARRAY_TASK_ID in $RUN_DIR"
python {project_root}/scr/simulation.py --run-dir "$RUN_DIR"
"""
        script_path = exp_dir / "submit.sh"
        with open(script_path, "w") as f:
            f.write(submit_script)

    elif runner == "local":
        workers = run_config.get("workers", 4)
        local_script = f"""#!/bin/bash
# Local parallel execution using xargs
echo "Running {total_jobs} jobs locally across {workers} workers..."

# Skip header, grab column 2 (run_dir), and pass to xargs
tail -n +2 {manifest_path} | cut -d',' -f2 | xargs -I {{}} -P {workers} bash -c '
    echo "Starting job in {{}}"
    python {project_root}/src/simulation.py --run-dir "$RUN_DIR"
'
"""
        script_path = exp_dir / "run_local.sh"
        with open(script_path, "w") as f:
            f.write(local_script)
        # Make the local script executable
        script_path.chmod(0o755)

    print(f"✅ Generated {total_jobs} jobs in {exp_dir}")
    print(f"📄 Manifest written to: {manifest_path}")
    print(f"🚀 Execution script: {script_path}")


if __name__ == "__main__":
    main()
