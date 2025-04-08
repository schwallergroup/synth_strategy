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

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis route uses a late-stage amide coupling strategy,
    where an amide bond is formed in one of the final steps of the synthesis.
    """
    amide_coupling_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_depth

        if node["type"] == "reaction":
            # Check if this is an amide coupling reaction
            if "rsmi" in node.get("metadata", {}):
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for amide coupling reaction types
                    is_amide_coupling = False

                    # Check for specific amide coupling reaction types
                    amide_reaction_types = [
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Carboxylic acid with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Schotten-Baumann_amide",
                    ]

                    for reaction_type in amide_reaction_types:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Detected {reaction_type} at depth {depth}")
                            is_amide_coupling = True
                            break

                    # If no specific reaction type matched, check for functional group changes
                    if not is_amide_coupling:
                        # Check for carboxylic acid in reactants
                        acid_present = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants
                        )

                        # Check for acyl halides (alternative to carboxylic acid)
                        acyl_halide_present = any(
                            checker.check_fg("Acyl halide", r) for r in reactants
                        )

                        # Check for amine in reactants
                        primary_amine_present = any(
                            checker.check_fg("Primary amine", r) for r in reactants
                        )
                        secondary_amine_present = any(
                            checker.check_fg("Secondary amine", r) for r in reactants
                        )

                        # Check for amide in product
                        primary_amide_in_product = checker.check_fg("Primary amide", product)
                        secondary_amide_in_product = checker.check_fg("Secondary amide", product)
                        tertiary_amide_in_product = checker.check_fg("Tertiary amide", product)

                        amine_present = primary_amine_present or secondary_amine_present
                        amide_in_product = (
                            primary_amide_in_product
                            or secondary_amide_in_product
                            or tertiary_amide_in_product
                        )

                        if (
                            (acid_present or acyl_halide_present)
                            and amine_present
                            and amide_in_product
                        ):
                            print(
                                f"Amide coupling detected through functional group analysis at depth {depth}"
                            )
                            print(f"  Acid present: {acid_present}")
                            print(f"  Acyl halide present: {acyl_halide_present}")
                            print(f"  Amine present: {amine_present}")
                            print(f"  Amide in product: {amide_in_product}")
                            is_amide_coupling = True

                    if is_amide_coupling:
                        amide_coupling_depth = depth
                        print(f"Amide coupling confirmed at depth {depth}")

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if amide coupling occurred in the first two steps (late stage)
    if amide_coupling_depth is not None and amide_coupling_depth <= 1:
        print(f"Late-stage amide coupling strategy detected at depth {amide_coupling_depth}")
        return True

    print("No late-stage amide coupling detected")
    return False
