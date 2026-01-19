import os

# analysis/model_select.py

def find_pdb_files(root):
    pdbs = []
    for dirpath, _, filenames in os.walk(root):
        for f in filenames:
            if f.endswith(".pdb"):
                pdbs.append(os.path.join(dirpath, f))
    return pdbs


def load_representative_model(tmpdir, tool):
    pdb_files = find_pdb_files(tmpdir)

    if not pdb_files:
        raise ValueError("No PDB files found in zip")

    if tool == "boltz":
        # e.g. model_0.pdb
        pdb_files.sort(key=lambda x: ("model_0" not in x, x))

    # ?? TODO: schrodinger uses other file dir. (x pdb -> maegz)
    elif tool == "schrodinger":
        # e.g. poseviewer, *_pv.pdb
        pdb_files.sort(key=lambda x: ("pv" not in x.lower(), x))

    else:
        raise ValueError(f"Unknown tool: {tool}")
    
    return pdb_files[0]