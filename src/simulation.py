import hoomd
import hoomd.md
import hoomd.write
import numpy as np
import yaml
from pathlib import Path

project_root = Path.cwd().parent
def load_hydropathy(scale_name):
    yaml_path = project_root / "configs" / "simulation" / "hydropathy_scales.yaml"
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    return data[scale_name]

AA_ORDER = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
hps_dict = load_hydropathy("HPS1")
missing = set(AA_ORDER) - set(hps_dict.keys())
if missing:
    raise ValueError(f"Hydropathy scale is missing values for: {missing}")
    
ids = {aa: i for i, aa in enumerate(AA_ORDER)}

hps = np.array([hps_dict[aa] for aa in AA_ORDER])

masses = np.array([71.08, 156.2, 114.1, 115.1, 103.1, 128.1, 129.1, 57.05, 137.1, 113.2
, 113.2, 128.2, 131.2, 147.2, 97.12, 87.08, 101.1, 186.2, 163.2, 99.07], dtype=np.float32)

charges = np.array([0.0, 1.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0
, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float32)

sigmas = np.array([0.504, 0.656, 0.568, 0.558, 0.548, 0.602, 0.592, 0.45, 0.608, 0.618, 0.618, 0.636
 , 0.618, 0.636, 0.556, 0.518, 0.562, 0.678, 0.646, 0.586], dtype=np.float32)

sequence  = 'GSHCFLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDC'
N         = len(sequence)
GSD_FILE  = 'very_cold_run.gsd'
BOX_SIZE  = 500    # nm, appropriate for a 23-residue chain
T_Kelvin  = 5      #  ~room temp in HPS
DT        = 0.01     # timestep
N_STEPS   = int((5*(N)**2.2)/DT)
SAVE_EVERY= 1000      # save frame every N steps -> 100 frames total
EPS_HPS   = 0.2       # kJ/mol, fixed LJ epsilon for all pairs
# ──────────────────────────────────────────────────────────────────────────────
KT = T_Kelvin*0.001987204259 # Boltzmann constant in kCal/mol/K
EPSILON = 0.1 # KCal/mol

sqid = np.array([ids[aa] for aa in sequence], dtype=np.int32)

seq_mass = masses[sqid]
seq_charge = charges[sqid]
seq_sigma = sigmas[sqid]
seq_hp = hps[sqid]

ionic_concentration = 150*1e-3 # in M or mol/L

fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3 #temperature dependent dielectric constant of water
epsw = fepsw(T_Kelvin) # dielectric constant of water at T 
lB = (1.6021766**2/(4*np.pi*8.854188*epsw))*(6.022*1000/KT)/4.184 # Bjerrum length in nm

yukawa_eps = lB*KT
yukawa_kappa = np.sqrt(8*np.pi*lB*ionic_concentration*6.022/10)

# ── Build snapshot ────────────────────────────────────────────────────────────
device = hoomd.device.auto_select()
sim = hoomd.Simulation(device=device, seed=42)

snapshot = hoomd.Snapshot()
snapshot.configuration.box = [BOX_SIZE, BOX_SIZE, BOX_SIZE, 0, 0, 0]
snapshot.particles.types = AA_ORDER
snapshot.particles.N = N
snapshot.bonds.N = N - 1      # must pre-allocate before populate_snapshot
snapshot.bonds.types = ['AA-bond']
# Particles
snapshot.particles.typeid[:] = sqid
snapshot.particles.mass[:]   = seq_mass
snapshot.particles.charge[:] = seq_charge
snapshot.particles.diameter[:] = seq_sigma

# Place beads in a straight line with ~0.38 nm spacing (realistic Cα-Cα)
snapshot.particles.position[:] = [[i * 0.38 - (N * 0.38) / 2, 0, 0] for i in range(N)]

# Bonds: chain connectivity
snapshot.bonds.N = N - 1
snapshot.bonds.types = ['AA-bond']
snapshot.bonds.typeid[:] = [0] * (N - 1)
snapshot.bonds.group[:]  = [[i, i + 1] for i in range(N - 1)]

sim.create_state_from_snapshot(snapshot)

# ── Forces ────────────────────────────────────────────────────────────────────
nl = hoomd.md.nlist.Tree(buffer=0.2)
nl.exclusions = ['bond']
lj1 = hoomd.md.pair.LJ(nlist=nl, mode='shift')
lj2 = hoomd.md.pair.LJ(nlist=nl)


yukawa = hoomd.md.pair.Yukawa(nlist=nl, default_r_cut=3.5, mode='shift')

for i, aa1 in enumerate(AA_ORDER):
    for j, aa2 in enumerate(AA_ORDER):
        if j < i:
            continue
        sig_ij = (sigmas[i] + sigmas[j]) / 2
        lam_ij = (hps[i]   + hps[j])   / 2
        lj2.params[(aa1, aa2)] = dict(epsilon=EPSILON*lam_ij, sigma=sig_ij)
        lj1.params[(aa1, aa2)] = dict(epsilon=EPSILON*(1-lam_ij), sigma=sig_ij)
        
        lj1.r_cut[(aa1, aa2)] = 2**(1/6) * sig_ij   # purely repulsive (WCA)

        lj2.r_cut[(aa1, aa2)] = 3.5        # attractive
        yukawa.params[(aa1, aa2)] = dict(epsilon=yukawa_eps*charges[i] * charges[j], kappa=yukawa_kappa)
        yukawa.r_cut[(aa1, aa2)]  = 3.5

harmonic = hoomd.md.bond.Harmonic()
harmonic.params['AA-bond'] = dict(k=1920, r0=0.38)
# ── Integrator ────────────────────────────────────────────────────────────────
langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=KT)
sim.operations.integrator = hoomd.md.Integrator(
    dt=DT,
    methods=[langevin],
    forces=[lj1, lj2, yukawa, harmonic]
)

# ── Writer ────────────────────────────────────────────────────────────────────
gsd_writer = hoomd.write.GSD(
    filename=GSD_FILE,
    trigger=hoomd.trigger.Periodic(SAVE_EVERY),
    mode='wb'
)
sim.operations.writers.append(gsd_writer)

print(f"Running HPS simulation for {N_STEPS} steps...")
print(f"Sequence: {sequence} ({N} residues)")
print(f"Saving every {SAVE_EVERY} steps -> {N_STEPS // SAVE_EVERY} frames to {GSD_FILE}\n")

sim.run(N_STEPS)
# Force flush to disk
gsd_writer.flush()
print("Wrote", GSD_FILE)
print("Simulation complete!")
