from rdkit import Chem
from rdkit.Chem import rdMolAlign


def _load_mol(path, sanitize=True, removeHs=False):
    mol = Chem.MolFromPDBFile(
        path,
        sanitize=sanitize,
        removeHs=removeHs
    )
    return mol


def compute_rdkit_ligand_rmsd(
    ref_lig_pdb: str,
    mob_lig_pdb: str,
) -> float:
    # --- Attempt 1: strict ---
    try:
        ref = _load_mol(ref_lig_pdb, sanitize=True, removeHs=False)
        mob = _load_mol(mob_lig_pdb, sanitize=True, removeHs=False)

        if ref is None or mob is None:
            raise ValueError("Mol load failed")

        match = mob.GetSubstructMatch(ref)
        if not match:
            raise ValueError("Substruct match failed")

        atom_map = [(i, j) for i, j in enumerate(match)]
        return rdMolAlign.CalcRMS(mob, ref, atomMap=atom_map)

    except Exception as e1:
        pass

    # --- Attempt 2: symmetry-aware RMSD ---
    try:
        return rdMolAlign.GetBestRMS(mob, ref)
    except Exception:
        pass

    # --- Attempt 3: no sanitize ---
    try:
        ref = _load_mol(ref_lig_pdb, sanitize=False, removeHs=False)
        mob = _load_mol(mob_lig_pdb, sanitize=False, removeHs=False)
        return rdMolAlign.GetBestRMS(mob, ref)
    except Exception:
        pass

    # --- Attempt 4: heavy atoms only ---
    try:
        ref = _load_mol(ref_lig_pdb, sanitize=False, removeHs=True)
        mob = _load_mol(mob_lig_pdb, sanitize=False, removeHs=True)
        return rdMolAlign.GetBestRMS(mob, ref)
    except Exception:
        pass

    raise RuntimeError(
        "RDKit RMSD failed after all fallbacks "
        "(likely different chemistry or binding mode)"
    )
