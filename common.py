"""Settings and ligand preparation shared by docking.py and validate_redock.py."""
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

IVOSIDENIB = "Fc1cncc(c1)N([C@H](C(=O)NC1CC(F)(F)C1)c1ccccc1Cl)C(=O)[C@@H]1CCC(=O)N1c1cc(ccn1)C#N"
VORASIDENIB = "C[C@@H](Nc1nc(N[C@H](C)C(F)(F)F)nc(n1)c2cccc(Cl)n2)C(F)(F)F"  # RCSB ligand 9UO
BOX_SIZE = [22.0, 22.0, 22.0]
EXHAUSTIVENESS = 16
SEEDS = [1, 2, 3]


def ligand_pdbqt(smiles):
    """SMILES -> 3D conformer (ETKDG + MMFF) -> PDBQT string via Meeko."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    setup = MoleculePreparation().prepare(mol)[0]
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise ValueError(err)
    return pdbqt
