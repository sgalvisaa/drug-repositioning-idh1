"""Ligand-based screen of the Broad Drug Repurposing Hub for IDH1 R132H glioma candidates.

Funnel: approved/launched drugs -> valid SMILES -> relaxed Lipinski -> BBB filter
-> top 10 by Morgan-fingerprint Tanimoto similarity to ivosidenib.

Usage: python screening.py   (expects the two Hub files in data/, see README)
"""
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

drug_info = pd.read_csv("data/repo-drug-annotation-20200324.txt", sep="\t", comment="!")
sample_info = pd.read_csv("data/repo-sample-annotation-20240610.txt", sep="\t", comment="!")
print(f"Hub drugs: {len(drug_info):,} | sample rows: {len(sample_info):,}")

merged = pd.merge(sample_info, drug_info, on="pert_iname", how="inner").dropna(subset=["smiles"])
approved = merged[merged["clinical_phase"].isin(["Phase 4", "Launched"])]
print(f"Approved/launched sample rows: {len(approved):,}")
approved = approved.drop_duplicates(subset="pert_iname").copy()
print(f"Unique approved drugs: {len(approved):,}")

approved["mol"] = approved["smiles"].apply(Chem.MolFromSmiles)
approved = approved.dropna(subset=["mol"]).copy()
print(f"Valid SMILES: {len(approved):,}")

approved["MW"] = approved["mol"].apply(Descriptors.MolWt)
approved["LogP"] = approved["mol"].apply(Descriptors.MolLogP)
approved["HBD"] = approved["mol"].apply(rdMolDescriptors.CalcNumHBD)
approved["HBA"] = approved["mol"].apply(rdMolDescriptors.CalcNumHBA)
approved["TPSA"] = approved["mol"].apply(Descriptors.TPSA)

# Relaxed Lipinski (MW <= 600) because ivosidenib itself is ~583 Da
lipinski = approved[(approved.MW <= 600) & (approved.LogP <= 5) & (approved.HBD <= 5)
                    & (approved.HBA <= 10) & (approved.TPSA <= 140)].copy()
print(f"Relaxed Lipinski: {len(lipinski):,}")

ivo_mol = approved.loc[approved.pert_iname == "ivosidenib", "mol"].iloc[0]
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
ivo_fp = fpgen.GetFingerprint(ivo_mol)
lipinski["tanimoto"] = lipinski["mol"].apply(
    lambda m: DataStructs.TanimotoSimilarity(ivo_fp, fpgen.GetFingerprint(m)))

# CNS penetration heuristic for glioma: TPSA < 90, MW < 500, LogP > 1
bbb = lipinski[(lipinski.TPSA < 90) & (lipinski.MW < 500) & (lipinski.LogP > 1)]
print(f"BBB filter: {len(bbb):,}")

candidates = (bbb[bbb.pert_iname != "ivosidenib"]
              .sort_values("tanimoto", ascending=False).head(10))
print(f"\nTop 10 by Tanimoto to ivosidenib (BBB-passing pool max = {bbb.tanimoto.max():.3f}):")
print(candidates[["pert_iname", "tanimoto", "MW", "LogP", "TPSA", "moa"]].round(3).to_string(index=False))

cols = ["pert_iname", "smiles", "tanimoto", "MW", "LogP", "TPSA", "moa"]
candidates[cols].round(3).to_csv("results/docking_candidates.csv", index=False)
# Random BBB-passing drugs (not in the top 10) as a docking background for comparison
bbb[~bbb.pert_iname.isin(candidates.pert_iname) & (bbb.pert_iname != "ivosidenib")] \
    .sample(40, random_state=7)[cols].round(3).to_csv("results/background_drugs.csv", index=False)
