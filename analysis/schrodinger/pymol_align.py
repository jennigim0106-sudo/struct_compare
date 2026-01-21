# analysis/schrodinger/pymol_align.py


def align_schrodinger_protein_to_reference(
    cmd,
    group_name: str,
    ref_obj: str = "exp",
) -> float:

    mob_sel = f"{group_name} and polymer.protein"
    ref_sel = f"{ref_obj} and polymer.protein"

    rmsd = cmd.align(mob_sel, ref_sel)[0]

    return rmsd