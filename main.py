import streamlit as st
import pymol2
import os
import tempfile
import pandas as pd

from analysis import pymol_ligand_rmsd as plr
from analysis.io import fetch_pdb_structure


RESULT_DIR = "./results/interactions"
os.makedirs(RESULT_DIR, exist_ok=True)



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
    "Upload Schrödinger glide-dock result (.zip)",
    type=["zip"],
    accept_multiple_files=False
)

# --- Schrödinger docking pose selection ---
schro_pose_selection = None

if schro_zip is not None:
    import tempfile
    from analysis.schrodinger_io import extract_schrodinger_zip, parse_glide_csv

    if "schro_tmpdir" not in st.session_state:
        st.session_state["schro_tmpdir"] = tempfile.mkdtemp()

    tmpdir = st.session_state["schro_tmpdir"]

    zip_path = os.path.join(tmpdir, schro_zip.name)
    with open(zip_path, "wb") as f:
        f.write(schro_zip.read())

    extract_dir = extract_schrodinger_zip(zip_path, tmpdir)
    df = parse_glide_csv(extract_dir)

    st.subheader("Schrödinger docking poses (click to select)")

    df_display = df[["pose_id", "glide_score"]].copy()
    df_display["select"] = False

    edited_df = st.data_editor(
        df_display,
        hide_index=True,
        column_config={
            "select": st.column_config.CheckboxColumn(
                "Select",
                help="Click to include/exclude this pose"
            )
        },
        disabled=["pose_id", "glide_score"],
        use_container_width=True
    )

    selected_poses = edited_df.loc[edited_df["select"], "pose_id"].tolist()


    schro_pose_selection = {
        "zip_path": zip_path,
        "extract_dir": extract_dir,
        "selected_poses": selected_poses
    }




# 3. Run Analysis button
if st.button("Run Analysis"):
    if not pdb_id:
        st.error("Please enter a PDB ID")
        st.stop()

    if (schro_zip is not None) and (not schro_pose_selection):
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
        if boltz_zips:

            boltz_tmpdir = tempfile.mkdtemp(prefix="boltz_")
            boltz_zip_paths = []

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
                tool="boltz",
                result_dir=RESULT_DIR
            )
        st.session_state["results"] = results


# Representing Results

results = st.session_state.get("results")

if results and "Boltz" in results:
    st.subheader("Boltz Interaction Results")

    for run in results["Boltz"]["per_run"]:
        csv_path = run.get("interaction_csv")
        run_name = run.get("run_name", "Boltz")

        st.markdown(f"### {run_name}")

        # Protein alignment info
        st.write(f"Protein Alignment RMSD: {run['protein_rmsd']:.3f} Å")
        st.write(f"Aligned atoms: {run['aligned_atoms']}")

        # Ligand RMSD
        st.write(f"Ligand RMSD: {run['ligand_rmsd']:.3f}")

        # Interaction CSV table
        if csv_path and os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            st.dataframe(df)

            with open(csv_path, "rb") as f:
                st.download_button(
                    label=f"Download {run_name} interactions",
                    data=f,
                    file_name=os.path.basename(csv_path),
                    mime="text/csv",
                )
        else:
            st.warning(f"Interaction CSV not found: {csv_path}")


        # Plot images
        plot_files = run.get("plot_files", [])
        plot_files = sorted(
            plot_files,
            key=lambda p: (
                "affinity" not in p,
                "pae" not in p,
                "plddt" not in p
            )
        )
        
        if plot_files:
            st.subheader("Associated Plots")

            # 1) 첫 번째 이미지 단독 표시 (affinity)
            st.image(plot_files[0], width=700)

            # 2) 나머지 이미지들 (PAE, pLDDT) → 한 줄 2개
            if len(plot_files) > 1:
                cols = st.columns(2)
                for col, plot_path in zip(cols, plot_files[1:3]):
                    col.image(plot_path, width=350)

