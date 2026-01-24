# analysis/pymol_ligand_rmsd.py
import tempfile
import os
from pymol import cmd
from .io import extract_zip
from .model_select import load_representative_model
from .pymol_align import align_proteins
import gemmi
import pathlib

from analysis.pymol_interactions import extract_interactions


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

'''
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

'''

import gemmi

def parse_nonpoly_comp_ids(cif_file_path):
    # 1. CIF 파일 읽기
    doc = gemmi.cif.read(cif_file_path)
    block = doc.sole_block() # 또는 doc.find_block('data_block_name')

    # 2. _pdbx_entity_nonpoly 테이블 가져오기
    table = block.find('_pdbx_entity_nonpoly.', ['entity_id', 'name', 'comp_id'])

    nonpoly_info = []
    
    # 3. 데이터 반복 처리
    for row in table:
        entity_id = row[0]
        name = row[1]
        comp_id = row[2]
        nonpoly_info.append({
            'entity_id': entity_id,
            'name': name, 
            'comp_id': comp_id,
        })

    return nonpoly_info



# -----------------------------
# Experimental ligand detection
# -----------------------------
def detect_experimental_ligand(exp_ori):
    # Step 1: Collect all nonpolymer compound IDs
    base, _ = os.path.splitext(exp_ori)
    cif_file = base + ".cif"
    get_comp_ids = parse_nonpoly_comp_ids(cif_file)
    comp_ids = [r["comp_id"] for r in get_comp_ids]
    
    # Step 2: Exclude known non-ligands
    ligands = [c for c in sorted(comp_ids) if c not in STANDARD_AMINO_ACIDS + DNA_RESNAMES + RNA_RESNAMES + PTM_RESNAMES + NOT_LIGANDS]

    if not ligands:
        return None

    # Step 3: Return first ligand (can be first alphabetically or first in list)
    return ligands[0]

def find_actual_pymol_ligand_resn(cmd, exp_obj, cif_comp_id):
    """
    CIF comp_id (e.g. A1AQQ) → PyMOL에서 실제 존재하는 resn 반환
    """
    # PDB resn은 최대 4자 → prefix match
    candidates = cmd.get_model(f"{exp_obj} and not polymer").atom
    resns = sorted(set(a.resn for a in candidates))

    for r in resns:
        if cif_comp_id.startswith(r):
            return r

    return None


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
def analyze_single_run(cmd, exp_ori, exp_obj, zip_file, tool, result_dir):
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
        print("Protein alignment RMSD:", align_result["rmsd"], "/ aligned atoms:", align_result["aligned_atoms"])

        ######### selecting ligand for comparison ###########
        ref_lig = detect_experimental_ligand(exp_ori)
        if ref_lig is None:
            raise ValueError("No experimental ligand detected")

        actual_resn = find_actual_pymol_ligand_resn(cmd, exp_obj, ref_lig)
        if actual_resn is None:
            raise ValueError(
                f"Experimental ligand '{ref_lig}' not found in PyMOL object (resn truncation issue)"
            )

        ref_lig_sel = f"{exp_obj} and resn {actual_resn}"


        if tool == "boltz":
            mob_lig_sel = f"{mob_obj} and resn LIG"

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


        ####### RDKit ligand rmsd calc #########
        from .rdkit_ligand_rmsd import compute_rdkit_ligand_rmsd

        lig_rmsd = compute_rdkit_ligand_rmsd(
            lig_pdbs["ref_lig_pdb"],
            lig_pdbs["mob_lig_pdb"]
        )

        # Ligand–protein interaction analysis
        interaction_csv = os.path.join(
            result_dir,
            f"{run_tag}_interactions.csv"
        )
        interaction_csv = str(interaction_csv)

        df_interactions = extract_interactions(cmd=cmd,
            protein_sel=f"{mob_obj} and polymer.protein",
            ligand_sel=mob_lig_sel,
            output_csv=interaction_csv,
        )

        import shutil  # 파일 복사 위해 필요

        # --- Copy plots from TMPDIR to permanent results folder ---
        from analysis.io import extract_plots_from_extracted_zip

        plot_files = extract_plots_from_extracted_zip(
            extracted_root=tmpdir,
            run_tag=run_tag,
            save_root="results/interactions/plots"
        )


        cmd.delete(mob_obj)
        return {
            "run_name": getattr(zip_file, "name", os.path.basename(zip_file)),
            "ligand_rmsd": lig_rmsd,
            "interaction_df": df_interactions,
            "interaction_csv": interaction_csv,
            "protein_rmsd": align_result["rmsd"],   
            "aligned_atoms": align_result["aligned_atoms"],
            "plot_files": plot_files
        }



# -----------------------------
# Multiple runs
# -----------------------------
def analyze_runs(cmd, exp_ori, exp_obj, zip_files, tool, result_dir):
    run_results = []
    for zip_file in zip_files:
        r = analyze_single_run(cmd, exp_ori, exp_obj, zip_file, tool, result_dir)
        run_results.append(r)
    
    # mean_rmsd
    lig_rmsds = [r["ligand_rmsd"] for r in run_results]
    mean_rmsd = sum(lig_rmsds) / len(lig_rmsds)

    return {"per_run": run_results, "mean_RMSD": mean_rmsd}
