from pathlib import Path

import hoomd
import hoomd.md
import hoomd.write
import numpy as np
import yaml

from constants import AA_ORDER, charges, ids, masses, sigma


def load_hydropathy(scale_name):
    project_root = Path.cwd().parent
    yaml_path = project_root / "configs" / "simulation" / "hydropathy_scales.yaml"
    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f)
    return data[scale_name]


hps_dict = load_hydropathy("HPS1")
missing = set(AA_ORDER) - set(hps_dict.keys())
if missing:
    raise ValueError(f"Hydropathy scale is missing values for: {missing}")

hps = np.array([hps_dict[aa] for aa in AA_ORDER])


def populate_snapshot(snapshot, sequence):

    sqid = np.array([ids[aa] for aa in sequence], dtype=np.int32)
    N = len(sequence)

    # Particles

    snapshot.particles.N = N
    snapshot.particles.typeid[:] = sqid
    snapshot.particles.mass[:] = masses[sqid]
    snapshot.particles.charge[:] = charges[sqid]
    snapshot.particles.diameter[:] = sigma[sqid]

    # Place beads in a straight line with ~0.38 nm spacing (realistic Cα-Cα)
    snapshot.particles.position[:] = [
        [i * 0.38 - (N * 0.38) / 2, 0, 0] for i in range(N)
    ]

    # Bonds: chain connectivity
    snapshot.bonds.N = N - 1
    snapshot.bonds.types = ["AA-bond"]
    snapshot.bonds.typeid[:] = [0] * (N - 1)
    snapshot.bonds.group[:] = [[i, i + 1] for i in range(N - 1)]

    return snapshot


# ── Settings ──────────────────────────────────────────────────────────────────
SEQUENCE = "NLYIQWLKDGGPSSGRPPPS"
GSD_FILE = "trajectory.gsd"
BOX_SIZE = 30.0  # nm, appropriate for a 23-residue chain
KT = 1.0  # reduced units, ~room temp in HPS
DT = 0.005  # timestep
N_STEPS = 100000
SAVE_EVERY = 1000  # save frame every N steps -> 100 frames total
EPS_HPS = 0.2  # kJ/mol, fixed LJ epsilon for all pairs
# ──────────────────────────────────────────────────────────────────────────────

# ── Build snapshot ────────────────────────────────────────────────────────────
N = len(SEQUENCE)
device = hoomd.device.auto_select()
sim = hoomd.Simulation(device=device, seed=42)

snapshot = hoomd.Snapshot()
snapshot.configuration.box = [BOX_SIZE, BOX_SIZE, BOX_SIZE, 0, 0, 0]
snapshot.particles.types = AA_ORDER
populate_snapshot(snapshot, SEQUENCE)

sim.create_state_from_snapshot(snapshot)

# ── Forces ────────────────────────────────────────────────────────────────────
nl = hoomd.md.nlist.Tree(buffer=0.2)
lj = hoomd.md.pair.LJ(nlist=nl)

# HPS pair interactions: WCA (repulsive) if mean hydropathy < 0.5, else full LJ
for i, aa1 in enumerate(AA_ORDER):
    for j, aa2 in enumerate(AA_ORDER):
        if j < i:
            continue
        sig_ij = (sigma[i] + sigma[j]) / 2
        lam_ij = (hps[i] + hps[j]) / 2
        lj.params[(aa1, aa2)] = dict(epsilon=EPS_HPS, sigma=sig_ij)
        if lam_ij < 0.5:
            lj.r_cut[(aa1, aa2)] = 2 ** (1 / 6) * sig_ij  # purely repulsive (WCA)
        else:
            lj.r_cut[(aa1, aa2)] = 2.5 * sig_ij  # attractive

# FENE bonds for chain connectivity
fene = hoomd.md.bond.FENEWCA()
fene.params["AA-bond"] = dict(k=30.0, r0=1.5, epsilon=1.0, sigma=0.47, delta=0.0)

# ── Integrator ────────────────────────────────────────────────────────────────
langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=KT)
sim.operations.integrator = hoomd.md.Integrator(
    dt=DT, methods=[langevin], forces=[lj, fene]
)

# ── Writer ────────────────────────────────────────────────────────────────────
gsd_writer = hoomd.write.GSD(
    filename=GSD_FILE, trigger=hoomd.trigger.Periodic(SAVE_EVERY), mode="wb"
)
sim.operations.writers.append(gsd_writer)

print(f"Running HPS simulation for {N_STEPS} steps...")
print(f"Sequence: {SEQUENCE} ({N} residues)")
print(
    f"Saving every {SAVE_EVERY} steps -> {N_STEPS // SAVE_EVERY} frames to {GSD_FILE}\n"
)

sim.run(N_STEPS)
# Force flush to disk
gsd_writer.flush()
print("Wrote", GSD_FILE)
print("Simulation complete!")
