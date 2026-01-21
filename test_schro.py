# test_schro.py
import os
import tempfile
import pymol2
from analysis.io import fetch_pdb_structure
from analysis.schrodinger.io import (
    extract_schrodinger_zip,
    find_glide_csv,
    convert_maegz_to_pdb,
)
from analysis.schrodinger.parse import parse_glide_csv
from analysis.schrodinger.poses import (
    find_pv_maegz,
    get_docking_pose_objects,
    map_selected_poses,
)
from analysis.schrodinger.pymol_align import align_schrodinger_protein_to_reference

# ---------- 설정 ----------
pdb_id = "6NUQ"  # 실험 구조 PDB
schro_zip_path = r"C:\Users\User01\Desktop\intern\Target\STAT3\00_STAT3-small molecule\6NUQ (SI109)\6nuq_schrodinger_run\glide-dock_SP_6nuq.zip"  # Schrödinger ZIP 경로
selected_pose_ids = [1, 3]  # 테스트용 선택된 pose 번호 (i_i_glide_lignum)

# ---------- 실험 구조 다운로드 ----------
exp_file = fetch_pdb_structure(pdb_id, save_dir="./data/pdb")
print(f"Experimental PDB saved to: {exp_file}")

# ---------- Schrödinger ZIP 압축 해제 ----------
schro_workdir = extract_schrodinger_zip(schro_zip_path)
print(f"Schrödinger workdir: {schro_workdir}")

# ---------- CSV 읽기 ----------
csv_path = find_glide_csv(schro_workdir)
schro_df = parse_glide_csv(csv_path)
print(f"Glide CSV loaded: {csv_path}")
print(schro_df[["i_i_glide_lignum", "r_i_glide_gscore"]])

# ---------- PyMOL headless session ----------
with pymol2.PyMOL() as pymol:
    cmd = pymol.cmd

    # 1. 실험 구조 로드
    cmd.load(exp_file, "exp")

    # 2. .maegz -> .pdb 변환 후 로드
    pv_maegz = find_pv_maegz(schro_workdir)
    schro_pdb = convert_maegz_to_pdb(pv_maegz)
    cmd.load(schro_pdb, "schro")
    print(f"Converted Schrödinger PDB: {schro_pdb}")

    # 3. ligand 객체 수집 및 선택
    ligand_objects = get_docking_pose_objects(cmd, "schro")
    selected_pose_map = map_selected_poses(ligand_objects, selected_pose_ids)
    print(f"Selected pose objects: {selected_pose_map}")

    # 4. protein alignment (Schrödinger -> experimental)
    protein_rmsd = align_schrodinger_protein_to_reference(cmd, group_name="schro", ref_obj="exp")
    print(f"Protein RMSD: {protein_rmsd:.3f} Å")
