#!/bin/bash
# Local parallel execution using xargs
echo "Running 1 jobs locally across 1 workers..."

# Skip header, grab column 2 (run_dir), and pass to xargs
tail -n +2 /home/andres/Research/machine-learning-IDPs/runs/20260320_165810_test_local_2/manifest.csv | cut -d',' -f2 | xargs -I {} -P 1 bash -c '
    echo "Starting job in {}"
    python /home/andres/Research/machine-learning-IDPs/src/simulation.py --run-dir "$RUN_DIR"
'
