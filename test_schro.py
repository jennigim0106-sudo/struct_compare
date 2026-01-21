# test_schro.py
import os
import tempfile
import pymol2
from analysis.io import fetch_pdb_structure
from analysis.schrodinger.io import extract_schrodinger_zip, find_glide_csv, find_schrodinger_complex_pdb
from analysis.schrodinger.parse import parse_glide_csv
from analysis.schrodinger.poses import find_pv_maegz, parse_group_name_from_log, get_docking_pose_objects, map_selected_poses
from analysis.schrodinger.pymol_align import align_schrodinger_protein_to_reference

# ---------- 설정 ----------
pdb_id = "6NUQ"  # 실험 구조 PDB
schro_zip_path = r"C:\Users\User01\Desktop\intern\Target\STAT3\00_STAT3-small molecule\6NUQ (SI109)\6nuq_schrodinger_run\glide-dock_SP_6nuq.zip"  # Schrödinger ZIP 경로
selected_pose_ids = [1, 3]  # 테스트용 선택된 pose 번호 (i_i_glide_lignum)

# ---------- 실험 구조 다운로드 ----------
cif_file = fetch_pdb_structure(pdb_id, save_dir="./data/pdb")
pdb_file = cif_file.replace(".cif", ".pdb")  # convert to PDB if fetch returns CIF
print(f"CIF saved to: {cif_file}")
print(f"PDB saved to: {pdb_file}")

# ---------- Schrödinger ZIP 압축 해제 ----------
schro_workdir = extract_schrodinger_zip(schro_zip_path)
print(f"Schrödinger workdir: {schro_workdir}")

# ---------- CSV 읽기 ----------
csv_path = find_glide_csv(schro_workdir)
schro_df = parse_glide_csv(csv_path)
print(f"Glide CSV loaded: {csv_path}")

# ---------- PyMOL headless session ----------
with pymol2.PyMOL() as pymol:
    cmd = pymol.cmd

    # 1. 실험 구조 로드
    cmd.load(pdb_file, "exp")

    # 2. _pv.maegz 및 로그 파일에서 group_name 추출
    pv_maegz = find_pv_maegz(schro_workdir)
    group_name = parse_group_name_from_log(schro_workdir)
    print(f"pv.maegz file: {pv_maegz}")
    print(f"group_name: {group_name}")

    # 3. pv.maegz 로드
    cmd.load(pv_maegz)

    # 4. ligand 객체 수집 및 선택
    ligand_objects = get_docking_pose_objects(cmd, group_name)
    selected_pose_map = map_selected_poses(ligand_objects, selected_pose_ids)
    print(f"Selected pose objects: {selected_pose_map}")

    # 5. protein alignment (Schrödinger -> experimental)
    protein_rmsd = align_schrodinger_protein_to_reference(
        cmd,
        group_name=group_name,
        ref_obj="exp"
    )
    print(f"Protein RMSD: {protein_rmsd:.3f} Å")
