# analysis/pymol_ligand_rmsd.py
import tempfile
import os
from pymol import cmd
from .io import extract_zip
from .model_select import load_representative_model
from .pymol_align import align_proteins


# Ligand exclusion lists
STANDARD_AMINO_ACIDS = [
    "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL"
]
DNA_RESNAMES = ["DA","DT","DG","DC","DI"]
RNA_RESNAMES = ["A","U","G","C"]
PTM_RESNAMES = ["PTR","SEP","TPO"]
NOT_LIGANDS = ["HOH","DOD","GOL","PEG"]


def parse_nonpoly_comp_ids(cif_file):
    comp_ids = []

    in_loop = False
    headers = []
    comp_id_idx = None

    current_row = []
    multiline = False
    multiline_value = ""

    with open(cif_file, "r") as f:
        for raw_line in f:
            line = raw_line.strip()

            # loop start
            if line == "loop_":
                in_loop = True
                headers = []
                comp_id_idx = None
                current_row = []
                continue

            if not in_loop:
                continue

            # get headers
            if line.startswith("_pdbx_entity_nonpoly."):
                headers.append(line)
                if line == "_pdbx_entity_nonpoly.comp_id":
                    comp_id_idx = len(headers) - 1
                continue

            # if header ended,  and comp_id doesn't exist, continue iteration
            if comp_id_idx is None:
                continue

            # multiline field start
            if line.startswith(";") and not multiline:
                multiline = True
                multiline_value = ""
                continue

            # inside multiline field
            if multiline:
                if line.startswith(";"):
                    multiline = False
                    current_row.append(multiline_value.strip())
                else:
                    multiline_value += line + " "
                continue

            # datas
            if line:
                tokens = line.split()
                for tok in tokens:
                    current_row.append(tok)

            # row 완성 여부 체크
            if len(current_row) == len(headers):
                comp_ids.append(current_row[comp_id_idx])
                current_row = []

            # 다음 loop 시작 시 종료
            if line == "loop_":
                break

    return comp_ids



# -----------------------------
# Experimental ligand detection
# -----------------------------
def detect_experimental_ligand(exp_ori):
    # Step 1: Collect all nonpolymer compound IDs
    base, _ = os.path.splitext(exp_ori)
    cif_file = base + ".cif"
    comp_ids = parse_nonpoly_comp_ids(cif_file)
    
    # Step 2: Exclude known non-ligands
    ligands = [c for c in sorted(comp_ids) if c not in STANDARD_AMINO_ACIDS + DNA_RESNAMES + RNA_RESNAMES + PTM_RESNAMES + NOT_LIGANDS]

    if not ligands:
        return None

    # Step 3: Return first ligand (can be first alphabetically or first in list)
    return ligands[0]



def export_aligned_ligands(cmd, ref_obj: str, mob_obj: str, ref_lig_sel: str, mob_lig_sel: str, out_dir: str, tag: str):
    # Exports aligned ligands as PDB files for RDKit RMSD.
    ref_lig_pdb = os.path.join(out_dir, f"{tag}_ref_lig.pdb")
    mob_lig_pdb = os.path.join(out_dir, f"{tag}_mob_lig.pdb")

    # Safety: selection existence check
    if cmd.count_atoms(ref_lig_sel) == 0:
        raise ValueError(f"Empty ref ligand selection: {ref_lig_sel}")
    if cmd.count_atoms(mob_lig_sel) == 0:
        raise ValueError(f"Empty mob ligand selection: {mob_lig_sel}")

    cmd.save(ref_lig_pdb, ref_lig_sel)
    cmd.save(mob_lig_pdb, mob_lig_sel)

    return {
        "ref_lig_pdb": ref_lig_pdb,
        "mob_lig_pdb": mob_lig_pdb,
    }



# -----------------------------
# Single run analysis
# -----------------------------
def analyze_single_run(cmd, exp_ori, exp_obj, zip_file, tool):
    with tempfile.TemporaryDirectory() as tmpdir:
        extract_zip(zip_file, tmpdir)
        model_file = load_representative_model(tmpdir, tool)

        mob_obj = f"model_{os.path.splitext(os.path.basename(zip_file))[0]}"
        cmd.load(model_file, mob_obj)

        ######## PROTEIN ALIGNMENT ##########
        align_result = align_proteins(
            cmd = cmd,
            ref_obj=exp_obj,
            mob_obj=mob_obj,
            ref_sel="polymer.protein",
            mob_sel="polymer.protein",
            method="align"   # or "super" later if needed
        )
        print("Protein alignment RMSD:", align_result["rmsd"])

        ######### selecting ligand for comparison ###########
        ref_lig = detect_experimental_ligand(exp_ori)
        if ref_lig is None:
            raise ValueError("No experimental ligand detected")

        ref_lig_sel = (
            f"{exp_obj} and hetatm and not polymer.protein "
            f"and not resn {'+'.join(NOT_LIGANDS)}"
        )
        mob_lig_sel = f"{mob_obj} and not polymer.protein and not resn {'+'.join(NOT_LIGANDS)}"

        out_dir = tmpdir
        run_tag = os.path.splitext(os.path.basename(
            getattr(zip_file, "name", zip_file)
        ))[0]

        lig_pdbs = export_aligned_ligands(
            cmd=cmd,
            ref_obj=exp_obj,
            mob_obj=mob_obj,
            ref_lig_sel=ref_lig_sel,
            mob_lig_sel=mob_lig_sel,
            out_dir=out_dir,
            tag=run_tag
        )

        from .rdkit_ligand_rmsd import compute_rdkit_ligand_rmsd

        lig_rmsd = compute_rdkit_ligand_rmsd(
            lig_pdbs["ref_lig_pdb"],
            lig_pdbs["mob_lig_pdb"]
        )

        cmd.delete(mob_obj)
        return {
            "run_name": getattr(zip_file, "name", os.path.basename(zip_file)),
            "ligand_rmsd": lig_rmsd
        }



# -----------------------------
# Multiple runs
# -----------------------------
def analyze_runs(cmd, exp_ori, exp_obj, zip_files, tool):
    run_results = []
    for zip_file in zip_files:
        r = analyze_single_run(cmd, exp_ori, exp_obj, zip_file, tool)
        run_results.append(r)
    
    # median_rmsd?
    return {"per_run": run_results}
