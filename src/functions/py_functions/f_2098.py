#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects if the synthesis involves sequential functionalization of an aromatic ring.
    """
    # Track functionalization steps
    functionalization_steps = []

    # List of aromatic functionalization reactions to check
    aromatic_functionalizations = [
        "Aromatic nitration with HNO3",
        "Aromatic nitration with NO3 salt",
        "Aromatic nitration with NO2 salt",
        "Aromatic nitration with alkyl NO2",
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Aromatic fluorination",
        "Friedel-Crafts acylation",
        "Friedel-Crafts alkylation",
        "Reduction of nitro groups to amines",
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Aromatic hydroxylation",
        "Aromatic sulfonyl chlorination",
    ]

    # List of aromatic functional groups to check
    aromatic_fgs = [
        "Nitro group",
        "Aniline",
        "Aromatic halide",
        "Phenol",
        "Sulfonamide",
        "Aromatic thiol",
        "Triflate",
        "Tosylate",
    ]

    def dfs_traverse(node):
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Extract depth safely
                depth_match = re.search(
                    r"Depth: (\d+)",
                    node.get("metadata", {}).get("ID", "Depth: -1") or "Depth: -1",
                )
                depth = int(depth_match.group(1)) if depth_match else -1

                # Check for aromatic functionalization reactions
                for reaction_type in aromatic_functionalizations:
                    if checker.check_reaction(reaction_type, rsmi):
                        functionalization_steps.append((reaction_type, depth))
                        print(f"Detected {reaction_type} at depth {depth}")

                # Check for functional group changes on aromatic rings
                product_fgs = []
                reactant_fgs = []

                for fg in aromatic_fgs:
                    # Check if product has this functional group
                    if checker.check_fg(fg, product_smiles):
                        product_fgs.append(fg)

                    # Check if any reactant has this functional group
                    for reactant in reactants_smiles:
                        if checker.check_fg(fg, reactant):
                            reactant_fgs.append(fg)

                # Identify new functional groups in product
                new_fgs = [fg for fg in product_fgs if fg not in reactant_fgs]
                removed_fgs = [fg for fg in reactant_fgs if fg not in product_fgs]

                if new_fgs:
                    for fg in new_fgs:
                        functionalization_steps.append((f"Addition of {fg}", depth))
                        print(f"Detected addition of {fg} at depth {depth}")

                if removed_fgs:
                    for fg in removed_fgs:
                        functionalization_steps.append((f"Removal of {fg}", depth))
                        print(f"Detected removal of {fg} at depth {depth}")

                # Check for specific transformations
                if "Nitro group" in reactant_fgs and "Aniline" in product_fgs:
                    functionalization_steps.append(("Nitro reduction", depth))
                    print(f"Detected nitro reduction at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth (higher depth = earlier in synthesis)
    functionalization_steps.sort(key=lambda x: x[1], reverse=True)

    # Check if we have at least 2 sequential aromatic functionalizations
    has_strategy = len(functionalization_steps) >= 2

    if has_strategy:
        step_types = [step[0] for step in functionalization_steps]
        print(f"Aromatic functionalization sequence: {' → '.join(step_types)}")

    return has_strategy
