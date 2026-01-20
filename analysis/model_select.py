import os

# analysis/model_select.py

def find_files_by_ext(root, extensions):
    files = []
    for dirpath, _, filenames in os.walk(root):
        for f in filenames:
            if any(f.endswith(ext) for ext in extensions):
                files.append(os.path.join(dirpath, f))
    return files


def load_representative_model(tmpdir, tool):
    if tool == "boltz":
        files = find_files_by_ext(tmpdir, [".pdb"])
        if not files:
            raise ValueError("No PDB files found for boltz")

        # model_0 우선
        files.sort(key=lambda x: ("model_0" not in os.path.basename(x), x))
        return files[0]

    elif tool == "schrodinger":
        files = find_files_by_ext(tmpdir, [".mae", ".maegz"])
        if not files:
            raise ValueError("No MAE/MAEGZ files found for schrodinger")

        # poseviewer 우선
        files.sort(key=lambda x: ("pv" not in os.path.basename(x).lower(), x))
        return files[0]

    else:
        raise ValueError(f"Unknown tool: {tool}")
