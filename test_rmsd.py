import pymol2

from analysis import pymol_ligand_rmsd as plr
from analysis.io import fetch_pdb_structure

pdb_id = "6NJS"
boltz_zips = ["/Users/gimjimin/Desktop/2026_intern/05_compare structure/Target/STAT3/00_STAT3-small molecule/6NJS (SD36)/boltz_run/6njs_2.zip"]

def get_experimental_structure(pdb_id):
    save_dir = "./data/pdb"
    pdb_path = fetch_pdb_structure(pdb_id, save_dir=save_dir)
    return pdb_path


with pymol2.PyMOL() as pymol:
    cmd = pymol.cmd

    print("=== PyMOL session started ===")

    # --- Load experimental structure ---
    exp_pdb = get_experimental_structure(pdb_id)
    cmd.load(exp_pdb, "exp")

    print(f"Loaded experimental structure: {exp_pdb}")
    print("Total atoms in exp:", cmd.count_atoms("exp"))

    # --- Inspect experimental ligands ---
    het_resns = sorted({a.resn for a in cmd.get_model("exp and hetatm").atom})
    print("HETATM resns in exp:", het_resns)

    # =========================
    # Run Boltz analysis
    # =========================
    results = {}

    if boltz_zips:
        print("\n=== Running Boltz RMSD analysis ===")
        results["Boltz"] = plr.analyze_runs(
            cmd=cmd,
            exp_ori=exp_pdb,
            exp_obj="exp",
            zip_files=boltz_zips,
            tool="boltz"
        )


    print("\n=== FINAL RESULTS ===")
    print(results)

    print("\n=== PyMOL session finished ===")