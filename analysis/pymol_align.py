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