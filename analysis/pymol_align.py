# analysis/pymol_align.py



def align_proteins(cmd, ref_obj, mob_obj, ref_sel="polymer.protein", mob_sel="polymer.protein", method="align"):
    ref = f"{ref_obj} and {ref_sel}"
    mob = f"{mob_obj} and {mob_sel}"

    if method == "align":
        result = cmd.align(mob, ref)
    elif method == "super":
        result = cmd.super(mob, ref)
    else:
        raise ValueError(f"Unknown alignment method: {method}")

    return {"rmsd": result[0], "aligned_atoms": result[1]}


from pymol import cmd


def align_pose_protein_to_reference(
    pose_protein: str,
    ref_protein: str,
    method: str = "align"
):
    # --- sanity check ---
    if cmd.count_atoms(pose_protein) == 0:
        raise ValueError(f"No atoms found in pose_protein selection: {pose_protein}")

    if cmd.count_atoms(ref_protein) == 0:
        raise ValueError(f"No atoms found in ref_protein selection: {ref_protein}")

    # --- alignment ---
    if method == "align":
        rmsd = cmd.align(pose_protein, ref_protein)[0]

    elif method == "super":
        rmsd = cmd.super(pose_protein, ref_protein)[0]

    elif method == "cealign":
        # cealign returns dict, RMSD key differs
        result = cmd.cealign(ref_protein, pose_protein)
        rmsd = result.get("RMSD", None)

    else:
        raise ValueError(
            f"Unknown alignment method '{method}'. "
            "Choose from 'align', 'super', or 'cealign'."
        )

    return rmsd