import streamlit as st
import pymol2
import os
import tempfile
import pandas as pd

from analysis import pymol_ligand_rmsd as plr
from analysis.io import fetch_pdb_structure

from analysis.schrodinger.io import (
    extract_schrodinger_zip,
    find_glide_csv,
)
from analysis.schrodinger.parse import parse_glide_csv
from analysis.schrodinger.poses import (
    find_pv_maegz,
    parse_group_name_from_log,
    get_docking_pose_objects,
    map_selected_poses,
)
from analysis.schrodinger.pymol_align import (
    align_schrodinger_protein_to_reference,
)



@st.cache_data
def get_experimental_structure(pdb_id):
    save_dir = "./data/pdb"
    pdb_dir = fetch_pdb_structure(pdb_id, save_dir=save_dir)
    return pdb_dir
    

#0. Title

st.title("Protein Structure Comparison")

# 1. Experimental structure

pdb_id = st.text_input("Enter PDB ID (experimental structure)",  
                       placeholder="e.g. 6TLC", 
                       help="Experimental structure will be downloaded from RCSB PDB")

# 2. Model uploads
boltz_zips = st.file_uploader(
    "Upload Boltz prediction results (.zip)",
    type=["zip"],
    accept_multiple_files=True
)

schro_zip = st.file_uploader(
    "Upload Schrödinger result (.zip)",
    type=["zip"]
)

# Preprocessing schrodinger file for selection
schro_workdir = None
schro_df = None

if schro_zip is not None and "schro_workdir" not in st.session_state:
    st.session_state["schro_workdir"] = extract_schrodinger_zip(schro_zip)
    csv_path = find_glide_csv(st.session_state["schro_workdir"])
    st.session_state["schro_df"] = parse_glide_csv(csv_path)

schro_workdir = st.session_state.get("schro_workdir")
schro_df = st.session_state.get("schro_df")

# pose selection
selected_pose_ids = []

if schro_df is not None:
    st.subheader("Schrödinger Docking Poses")

    st.dataframe(
        schro_df[["i_i_glide_lignum", "r_i_glide_gscore"]]
    )

    selected_pose_ids = st.multiselect(
        "Select docking poses",
        options = schro_df["i_i_glide_lignum"].astype(int).tolist()
    )



# 3. Run Analysis button
if st.button("Run Analysis"):
    if not pdb_id:
        st.error("Please enter a PDB ID")
        st.stop()


    if (schro_zip is not None) and (not selected_pose_ids):
        st.error("Please select at least one docking pose")
        st.stop()


    # PyMOL headless session
    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd

        # Experimental structure load
        exp_file = get_experimental_structure(pdb_id)
        cmd.load(exp_file, "exp")

        results = {}

        # Boltz runs
        boltz_zip_paths = []

        if boltz_zips:
            boltz_tmpdir = tempfile.mkdtemp(prefix="boltz_")

            for uf in boltz_zips:
                zip_path = os.path.join(boltz_tmpdir, uf.name)
                with open(zip_path, "wb") as f:
                    f.write(uf.read())
                boltz_zip_paths.append(zip_path)

            results["Boltz"] = plr.analyze_runs(
                cmd,
                exp_file,
                "exp",
                boltz_zip_paths,
                tool="boltz"
            )


        # Schrodinger runs

        if schro_workdir is not None:
            pv_maegz = find_pv_maegz(schro_workdir)
            group_name = parse_group_name_from_log(schro_workdir)

            # load pv structure
            cmd.load(pv_maegz)

            # collect ligand pose objects
            ligand_objects = get_docking_pose_objects(cmd, group_name)

            # map selected poses
            selected_pose_map = map_selected_poses(
                ligand_objects,
                selected_pose_ids
            )

            # initial settings for results dictionary
            results["Schrodinger"] = {
                "pv_maegz": pv_maegz,
                "group_name": group_name,
                "selected_poses": selected_pose_map,
            }

            # protein alignment (Schrödinger -> experimental)
            protein_rmsd = align_schrodinger_protein_to_reference(
                cmd,
                group_name=group_name,
                ref_obj="exp",
            )
            results["Schrodinger"]["protein_rmsd"] = protein_rmsd

    
    st.write(results)