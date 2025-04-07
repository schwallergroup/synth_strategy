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
    This function detects if the synthesis involves formation of a urea moiety
    in the final or penultimate step.
    """
    urea_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal urea_formation_detected

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            print(f"Checking reaction step at depth {depth}")
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check if any reactant contains urea or thiourea
                    reactants_list = reactants_smiles.split(".")
                    reactants_have_urea = any(
                        checker.check_fg("Urea", reactant)
                        for reactant in reactants_list
                    )
                    reactants_have_thiourea = any(
                        checker.check_fg("Thiourea", reactant)
                        for reactant in reactants_list
                    )

                    # Check if product contains urea or thiourea
                    product_has_urea = checker.check_fg("Urea", product_smiles)
                    product_has_thiourea = checker.check_fg("Thiourea", product_smiles)

                    print(f"Depth {depth} - Reactants have urea: {reactants_have_urea}")
                    print(
                        f"Depth {depth} - Reactants have thiourea: {reactants_have_thiourea}"
                    )
                    print(f"Depth {depth} - Product has urea: {product_has_urea}")
                    print(
                        f"Depth {depth} - Product has thiourea: {product_has_thiourea}"
                    )

                    # Check if this is a urea/thiourea formation reaction
                    is_urea_synthesis = (
                        checker.check_reaction(
                            "Urea synthesis via isocyanate and primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Urea synthesis via isocyanate and secondary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Urea synthesis via isocyanate and diazo", rsmi
                        )
                        or checker.check_reaction(
                            "Urea synthesis via isocyanate and sulfonamide", rsmi
                        )
                        or checker.check_reaction("urea", rsmi)
                        or checker.check_reaction("thiourea", rsmi)
                    )

                    print(
                        f"Depth {depth} - Is urea/thiourea synthesis reaction: {is_urea_synthesis}"
                    )

                    # Urea/thiourea formation is detected if:
                    # 1. Product has urea/thiourea but reactants don't, OR
                    # 2. The reaction is specifically a urea/thiourea synthesis
                    if (
                        (product_has_urea and not reactants_have_urea)
                        or (product_has_thiourea and not reactants_have_thiourea)
                        or is_urea_synthesis
                    ):
                        print(f"Urea/thiourea formation detected at depth {depth}")
                        urea_formation_detected = True
                except Exception as e:
                    print(
                        f"Error processing SMILES in urea_formation_strategy: {str(e)}"
                    )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: urea_formation_detected = {urea_formation_detected}")
    return urea_formation_detected
