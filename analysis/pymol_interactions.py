import pandas as pd
import os


DEFAULT_CUTOFFS = {
    "hbond": 3.5,
    "salt_bridge": 4.0,
    "pi_pi": 5.0,
    "hydrophobic": 4.0,
}


def _pair_to_record(cmd, prot_pair, lig_pair, interaction_type):
    (prot_obj, prot_idx) = prot_pair
    (lig_obj, lig_idx) = lig_pair

    prot_model = cmd.get_model(f"{prot_obj} and index {prot_idx}")
    lig_model  = cmd.get_model(f"{lig_obj} and index {lig_idx}")

    if len(prot_model.atom) == 0 or len(lig_model.atom) == 0:
        return None

    p = prot_model.atom[0]
    l = lig_model.atom[0]

    distance = cmd.get_distance(
        f"{prot_obj} and index {prot_idx}",
        f"{lig_obj} and index {lig_idx}"
    )

    return {
        "interaction_type": interaction_type,
        "protein_resn": p.resn,
        "protein_resi": p.resi,
        "protein_atom": p.name,
        "ligand_atom": l.name,
        "distance": round(distance, 3),
    }



def extract_interactions(cmd,
    protein_sel: str,
    ligand_sel: str,
    output_csv: str,
    cutoffs: dict = None,
):
    
    if cutoffs is None:
        cutoffs = DEFAULT_CUTOFFS

    records = []

    # protein_sel : e.g. "6njs and polymer.protein"
    # ligand_sel : e.g. mob_lig_sel = "6njs and resn LIG"

    # 1. Hydrogen bonds
    hbond_pairs = cmd.find_pairs(
        f"{protein_sel} and (donor or acceptor)",
        f"{ligand_sel} and (donor or acceptor)",
        cutoff=cutoffs["hbond"],
        mode=1,
    )

    for prot_atom, lig_atom in hbond_pairs:
        records.append(
            _pair_to_record(cmd, prot_atom, lig_atom, "hbond")
        )


    # 2. Salt bridges
    protein_pos = f"{protein_sel} and resn ARG+LYS+HIS"
    protein_neg = f"{protein_sel} and resn ASP+GLU"
    ligand_pos = f"{ligand_sel} and (formal_charge > 0)"
    ligand_neg = f"{ligand_sel} and (formal_charge < 0)"

    salt_pairs_1 = cmd.find_pairs(
        protein_pos,
        ligand_neg,
        cutoff=cutoffs["salt_bridge"],
        mode=1,
    )

    salt_pairs_2 = cmd.find_pairs(
        protein_neg,
        ligand_pos,
        cutoff=cutoffs["salt_bridge"],
        mode=1,
    )

    for prot_atom, lig_atom in salt_pairs_1 + salt_pairs_2:
        records.append(
            _pair_to_record(cmd, prot_atom, lig_atom, "salt_bridge")
        )


    # 3. π–π stacking (proxy)
    protein_aromatic = f"{protein_sel} and resn PHE+TYR+TRP+HIS"
    ligand_aromatic = f"{ligand_sel} and elem C"

    pi_pairs = cmd.find_pairs(
        protein_aromatic,
        ligand_aromatic,
        cutoff=cutoffs["pi_pi"],
        mode=1,
    )

    for prot_atom, lig_atom in pi_pairs:
        records.append(
            _pair_to_record(cmd, prot_atom, lig_atom, "pi_pi")
        )


    # 4. Hydrophobic contacts
    protein_hphob = (
        f"{protein_sel} and not (donor or acceptor) "
        "and not resn ASP+GLU+ARG+LYS+HIS"
    )
    
    ligand_hphob = (
        f"{ligand_sel} and not (donor or acceptor)"
    )

    hydrophobic_pairs = cmd.find_pairs(
        protein_hphob,
        ligand_hphob,
        cutoff=cutoffs["hydrophobic"],
        mode=1,
    )

    for prot_atom, lig_atom in hydrophobic_pairs:
        records.append(
            _pair_to_record(cmd, prot_atom, lig_atom, "hydrophobic")
        )

    # to csv
    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)

    return df