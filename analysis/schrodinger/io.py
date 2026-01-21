from pathlib import Path
import zipfile
import tempfile
import os

def extract_schrodinger_zip(zip_file_path: str) -> str:
    '''
    tmpdir = tempfile.mkdtemp(prefix="schro_")
    zip_path = os.path.join(tmpdir, zip_file.name)

    with open(zip_path,"wb") as f:
        f.write(zip_file.read())

    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(tmpdir)

    return tmpdir
    '''
    if not os.path.isfile(zip_file_path):
        raise FileNotFoundError(f"ZIP file not found: {zip_file_path}")

    tmpdir = tempfile.mkdtemp(prefix="schro_")

    with zipfile.ZipFile(zip_file_path, "r") as z:
        z.extractall(tmpdir)

    return tmpdir


def find_glide_csv(workdir: str) -> str:
    workdir = Path(workdir)

    csv_files = [
        p for p in workdir.rglob("glide-dock_*.csv")
        if not p.name.endswith("_skip.csv")
    ]

    if len(csv_files) == 0:
        raise RuntimeError("No suitable Glide CSV found")

    if len(csv_files) > 1:
        raise RuntimeError(
            f"Multiple candidate Glide CSVs found: "
            f"{[p.name for p in csv_files]}"
        )

    return str(csv_files[0])




# NOTE:
# This function is currently unused.
# Kept as a fallback utility in case PDB-based workflows are needed later.

def find_schrodinger_complex_pdb(workdir: str) -> str:
    """
    Find a Schrödinger-exported complex PDB inside extracted directory.
    Temporary rule:
      - first *.pdb file found
    """
    pdb_candidates = []

    for root, _, files in os.walk(workdir):
        for f in files:
            if f.lower().endswith(".pdb"):
                pdb_candidates.append(os.path.join(root, f))

    if not pdb_candidates:
        raise RuntimeError(
            "No PDB file found in Schrödinger result. "
            "Please export complex PDB from Maestro first."
        )

    return pdb_candidates[0]

