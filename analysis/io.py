# analysis/io.py
import requests
from Bio.PDB import MMCIFParser, PDBIO
import zipfile
import os

def fetch_pdb_structure(pdb_id, save_dir="data/pdb"):
    pdb_id = pdb_id.lower()
    os.makedirs(save_dir, exist_ok=True)

    cif_path = os.path.join(save_dir, f"{pdb_id}.cif")
    pdb_path = os.path.join(save_dir, f"{pdb_id}.pdb")

    # Download CIF if not already exists
    if not os.path.exists(cif_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        r = requests.get(url)
        if r.status_code != 200:
            raise ValueError(f"Failed to fetch PDB ID: {pdb_id}")
        with open(cif_path, "w") as f:
            f.write(r.text)
        print(f"CIF saved to: {cif_path}")
    else:
        print(f"CIF already exists: {cif_path}")

    # Convert CIF → PDB
    if not os.path.exists(pdb_path):
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(pdb_id, cif_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)
        print(f"PDB saved to: {pdb_path}")
    else:
        print(f"PDB already exists: {pdb_path}")

    return pdb_path


def extract_zip(zip_file, out_dir):
    with zipfile.ZipFile(zip_file, "r") as z:
        z.extractall(out_dir)
