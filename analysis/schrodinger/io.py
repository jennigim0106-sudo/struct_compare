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

import gzip
import shutil
import subprocess

def convert_maegz_to_pdb(maegz_path: str) -> str:

    mae_path = maegz_path.replace("_pv.maegz", "_pv.mae")
    pdb_path = mae_path.replace(".mae", ".pdb")

    # 1. decompress maegz -> mae
    with gzip.open(maegz_path, "rb") as f_in:
        with open(mae_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    # 2. mae -> pdb using structconvert
    subprocess.run(["structconvert", "-imae", mae_path, "-opdb", pdb_path], check=True)

    return pdb_path