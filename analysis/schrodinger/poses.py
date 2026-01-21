# analysis/schrodinger/poses.py

from pathlib import Path
import re



def find_pv_maegz(workdir: str) -> str:
    workdir = Path(workdir)
    maegz_files = list(workdir.rglob("*_pv.maegz"))

    if len(maegz_files) == 0:
        raise RuntimeError("No *_pv.maegz file found in Schrödinger results")

    if len(maegz_files) > 1:
        raise RuntimeError(
            f"Multiple *_pv.maegz files found: {[p.name for p in maegz_files]}"
        )

    return str(maegz_files[0])

################################

def parse_group_name_from_log(workdir: str) -> str:
    workdir = Path(workdir)
    log_files = list(workdir.rglob("*.log"))

    if not log_files:
        raise RuntimeError("No .log file found in Schrödinger results")

    pattern = re.compile(r"Writing\s+\d+\s+poses\s+to\s+(\S+)_pv\.maegz")

    for log_path in log_files:
        with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                m = pattern.search(line)
                if m:
                    return f"{m.group(1)}_pv"

    raise RuntimeError("Failed to parse group_name from Glide log file")


def get_docking_pose_objects(cmd, group_name: str) -> list[str]:
    members = cmd.get_object_list(group_name)

    ligand_objs = [
        obj for obj in members
        if cmd.count_atoms(f"{obj} and not polymer") > 0
    ]

    return ligand_objs


def map_selected_poses(
    ligand_objects: list[str],
    selected_pose_ids: list[int],
) -> dict[int, str]:
    pose_map = {}

    for lignum in selected_pose_ids:
        idx = lignum - 1  # Glide lignum is 1-based
        if idx < 0 or idx >= len(ligand_objects):
            raise IndexError(
                f"Selected lignum {lignum} out of range "
                f"(available poses: {len(ligand_objects)})"
            )
        pose_map[lignum] = ligand_objects[idx]

    return pose_map
