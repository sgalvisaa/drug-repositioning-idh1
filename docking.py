"""Structure-based docking of screening candidates into the IDH1 R132H allosteric pocket.

Receptor: PDB 8T7O (IDH1 R132H homodimer co-crystallized with ivosidenib), chains A + B.
Box: centered on the crystal ivosidenib, which sits at the dimer interface.

1. Positive controls: ivosidenib and vorasidenib (both approved mutant-IDH inhibitors).
2. The 10 screening candidates, 3 independent seeds each.
3. Background: 40 random drugs that passed the same filters, to show what a "typical" score is.

Pose validation against the crystal structure is in validate_redock.py.

Usage: python docking.py   (run screening.py first)
"""
import os
import subprocess

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from vina import Vina

from common import IVOSIDENIB, VORASIDENIB, BOX_SIZE, EXHAUSTIVENESS, SEEDS, ligand_pdbqt

RDLogger.DisableLog("rdApp.*")
os.makedirs("results/poses", exist_ok=True)

PDB = "structures/8T7O.pdb"

# --- Receptor: protein atoms of both chains (no waters, NADP, or ligand) ---
with open(PDB) as f:
    lines = f.readlines()
with open("structures/receptor_AB.pdb", "w") as f:
    f.writelines(l for l in lines if l.startswith("ATOM") and l[21] in "AB")
    f.write("END\n")
subprocess.run(["obabel", "structures/receptor_AB.pdb", "-xr", "-p", "7.4",
                "-O", "structures/receptor_AB.pdbqt"], check=True, capture_output=True)

# --- Crystal ivosidenib: coordinates from the PDB, bond orders from SMILES ---
with open("structures/ivosidenib_crystal.pdb", "w") as f:
    f.writelines(l for l in lines if l.startswith("HETATM") and l[17:20] == "IV3")
    f.write("END\n")
crystal = Chem.MolFromPDBFile("structures/ivosidenib_crystal.pdb")
crystal = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(IVOSIDENIB), crystal)
center = crystal.GetConformer().GetPositions().mean(axis=0).tolist()
print(f"Box center (crystal ivosidenib centroid): {np.round(center, 3)}, size {BOX_SIZE}")

# One Vina engine per seed; the affinity maps are computed once and reused for every ligand
engines = {}
for seed in SEEDS:
    engines[seed] = Vina(sf_name="vina", cpu=0, seed=seed, verbosity=0)
    engines[seed].set_receptor("structures/receptor_AB.pdbqt")
    engines[seed].compute_vina_maps(center=center, box_size=BOX_SIZE)


def run(name, smiles, seeds, group):
    rows = []
    heavy = Chem.MolFromSmiles(smiles).GetNumHeavyAtoms()
    pdbqt = ligand_pdbqt(smiles)
    for seed in seeds:
        v = engines[seed]
        v.set_ligand_from_string(pdbqt)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=9)
        if group != "background":
            v.write_poses(f"results/poses/{name}_seed{seed}.pdbqt", n_poses=9, overwrite=True)
        score = v.energies(n_poses=1)[0][0]
        rows.append(dict(compound=name, group=group, seed=seed, score=score,
                         heavy_atoms=heavy, ligand_efficiency=-score / heavy))
        print(f"  {group:<10} {name:<16} seed {seed}: {score:7.3f} kcal/mol")
    return rows


rows = []

# 1. Positive controls
print("\n[1] Ivosidenib (co-crystallized in 8T7O)")
rows += run("ivosidenib", IVOSIDENIB, SEEDS, "control")

print("\n[1] Vorasidenib (approved 2024 for IDH-mutant glioma)")
rows += run("vorasidenib", VORASIDENIB, SEEDS, "control")

# 2. Candidates
print("\n[2] Screening candidates")
for _, r in pd.read_csv("results/docking_candidates.csv").iterrows():
    rows += run(r.pert_iname, r.smiles, SEEDS, "candidate")

# 3. Background
print("\n[3] Background: random filter-passing drugs (seed 1 only)")
for _, r in pd.read_csv("results/background_drugs.csv").iterrows():
    try:
        rows += run(r.pert_iname, r.smiles, [1], "background")
    except Exception as e:  # e.g. embedding failure for an odd structure
        print(f"  skipped {r.pert_iname}: {e}")

df = pd.DataFrame(rows)
df.round(3).to_csv("results/docking_scores_all_seeds.csv", index=False)
summary = (df.groupby(["compound", "group"])
           .agg(best=("score", "min"), mean=("score", "mean"), sd=("score", "std"),
                heavy_atoms=("heavy_atoms", "first"), ligand_efficiency=("ligand_efficiency", "mean"))
           .reset_index().sort_values("mean"))
summary.round(3).to_csv("results/docking_summary.csv", index=False)
print("\n", summary[summary.group != "background"].round(2).to_string(index=False))
bg = summary[summary.group == "background"]["mean"]
print(f"\nBackground (n={len(bg)}): median {bg.median():.2f}, 90th pct {bg.quantile(0.1):.2f} kcal/mol")
