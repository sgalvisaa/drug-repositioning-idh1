"""Redock validation: can the docking setup put ivosidenib back where the crystal structure has it?

Two tests, same receptor, box, and Vina settings as docking.py:
  A. Self-docking: start from the crystal ligand's own geometry. Vina randomizes position,
     orientation, and rotatable torsions, but rings are rigid in Vina, so ring shapes are kept.
  B. From SMILES: start from an RDKit conformer, exactly how the candidates were prepared.
Pass = top-ranked pose within 2.0 A heavy-atom RMSD of the crystal pose (symmetry-aware, no re-alignment).

Usage: python validate_redock.py   (run docking.py first; it writes the receptor files)
"""
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
from vina import Vina

from common import IVOSIDENIB, BOX_SIZE, EXHAUSTIVENESS, SEEDS, ligand_pdbqt

RDLogger.DisableLog("rdApp.*")

crystal = AllChem.AssignBondOrdersFromTemplate(
    Chem.MolFromSmiles(IVOSIDENIB), Chem.MolFromPDBFile("structures/ivosidenib_crystal.pdb"))
center = crystal.GetConformer().GetPositions().mean(axis=0).tolist()
setup = MoleculePreparation().prepare(Chem.AddHs(crystal, addCoords=True))[0]
starts = {"A_self_docking_crystal_rings": PDBQTWriterLegacy.write_string(setup)[0],
          "B_from_smiles_rdkit_conformer": ligand_pdbqt(IVOSIDENIB)}

rows = []
for seed in SEEDS:
    v = Vina(sf_name="vina", cpu=0, seed=seed, verbosity=0)
    v.set_receptor("structures/receptor_AB.pdbqt")
    v.compute_vina_maps(center=center, box_size=BOX_SIZE)
    for test, pdbqt in starts.items():
        v.set_ligand_from_string(pdbqt)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=9)
        energies = v.energies(n_poses=9)[:, 0]
        poses = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule(v.poses(n_poses=9), skip_typing=True))[0])
        for conf in poses.GetConformers():
            rows.append(dict(test=test, seed=seed, pose_rank=conf.GetId() + 1,
                             score=energies[conf.GetId()],
                             rmsd=rdMolAlign.CalcRMS(poses, crystal, prbId=conf.GetId())))

df = pd.DataFrame(rows)
df.round(3).to_csv("results/redock_rmsd.csv", index=False)
top = df[df.pose_rank == 1].groupby("test").agg(score=("score", "mean"), rmsd_mean=("rmsd", "mean"),
                                                 rmsd_max=("rmsd", "max"))
top["pass_2A"] = top.rmsd_max < 2.0
print(top.round(2).to_string())

# Is the crystal pose itself the best-scoring arrangement? (scoring vs. sampling failure)
v.set_ligand_from_string(starts["A_self_docking_crystal_rings"])
print(f"\nCrystal pose: score_only {v.score()[0]:.2f}, after local optimization {v.optimize()[0]:.2f} kcal/mol")
