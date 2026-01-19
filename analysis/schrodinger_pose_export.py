import os
import subprocess

def export_selected_poses_to_pdb(
    pv_maegz_path: str,
    pose_ids: list[int],
    out_dir: str
):
    """
    Export selected poses from _pv.maegz to individual PDB files
    Requires Schrödinger structconvert in PATH
    """

    os.makedirs(out_dir, exist_ok=True)

    pdb_files = []

    for pid in pose_ids:
        out_pdb = os.path.join(out_dir, f"{pid}.pdb")

        cmd = [
            "structconvert",
            "-imae", pv_maegz_path,
            "-omae", out_pdb,
            "-n", str(pid)
        ]

        subprocess.run(cmd, check=True)
        pdb_files.append(out_pdb)

    return pdb_files
