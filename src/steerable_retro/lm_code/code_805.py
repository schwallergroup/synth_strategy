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
    This function detects if morpholine is introduced in the late stage of synthesis (depth 0-1)
    via nucleophilic substitution on an alkyl chain.
    """
    morpholine_introduced = False

    def dfs_traverse(node, current_depth=0):
        nonlocal morpholine_introduced

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Use current_depth from DFS traversal instead of trying to parse from ID
            if current_depth <= 1:  # Late stage (depth 0 or 1)
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Analyzing reaction at depth {current_depth}: {rsmi}")

                    # Check if morpholine is in the product
                    product_has_morpholine = checker.check_ring("morpholine", product)
                    if product_has_morpholine:
                        print(f"Product contains morpholine: {product}")

                    # Find which reactant is morpholine
                    morpholine_reactant = None
                    other_reactants = []

                    for reactant in reactants:
                        if checker.check_ring("morpholine", reactant):
                            morpholine_reactant = reactant
                            print(f"Found morpholine in reactant: {reactant}")
                        else:
                            other_reactants.append(reactant)

                    # Check if this is a nucleophilic substitution reaction
                    is_nucleophilic_sub = any(
                        [
                            checker.check_reaction(
                                "N-alkylation of secondary amines with alkyl halides", rsmi
                            ),
                            checker.check_reaction(
                                "N-alkylation of primary amines with alkyl halides", rsmi
                            ),
                            checker.check_reaction("Williamson Ether Synthesis", rsmi),
                            checker.check_reaction("Mitsunobu aryl ether", rsmi),
                            checker.check_reaction("S-alkylation of thiols", rsmi),
                            checker.check_reaction("S-alkylation of thiols with alcohols", rsmi),
                        ]
                    )

                    # If not identified by reaction type but morpholine appears in both reactant and product,
                    # it's likely a nucleophilic substitution
                    if not is_nucleophilic_sub and morpholine_reactant and product_has_morpholine:
                        print(
                            "Reaction not identified by type but morpholine appears in both reactant and product"
                        )
                        is_nucleophilic_sub = True

                    if is_nucleophilic_sub:
                        print(f"Found nucleophilic substitution reaction: {rsmi}")

                        # If morpholine is a reactant and product contains morpholine
                        if morpholine_reactant and product_has_morpholine:
                            # Check if the other reactant has a leaving group
                            for other_reactant in other_reactants:
                                has_leaving_group = any(
                                    [
                                        checker.check_fg("Primary halide", other_reactant),
                                        checker.check_fg("Secondary halide", other_reactant),
                                        checker.check_fg("Tertiary halide", other_reactant),
                                        checker.check_fg("Tosylate", other_reactant),
                                        checker.check_fg("Mesylate", other_reactant),
                                        checker.check_fg("Triflate", other_reactant),
                                        "Cl" in other_reactant,  # Check for chlorine atom
                                        "Br" in other_reactant,  # Check for bromine atom
                                        "I" in other_reactant,  # Check for iodine atom
                                    ]
                                )

                                if has_leaving_group:
                                    print(
                                        f"Found late-stage morpholine introduction at depth {current_depth}"
                                    )
                                    print(f"Reaction SMILES: {rsmi}")
                                    print(f"Morpholine reactant: {morpholine_reactant}")
                                    print(f"Other reactant with leaving group: {other_reactant}")
                                    morpholine_introduced = True
                                    break
                                else:
                                    print(
                                        f"No suitable leaving group found in reactant: {other_reactant}"
                                    )
                        elif morpholine_reactant:
                            print(f"Morpholine found in reactant but not in product: {product}")
                        elif product_has_morpholine:
                            print(f"Morpholine found in product but not in reactants")

                            # Check if morpholine is formed in the reaction (e.g., ring closure)
                            # This is a different mechanism than nucleophilic substitution
                            for reactant in reactants:
                                # Check if reactant contains components that could form morpholine
                                if checker.check_fg("Primary amine", reactant) and checker.check_fg(
                                    "Primary alcohol", reactant
                                ):
                                    print(
                                        f"Possible morpholine formation from components in reactant: {reactant}"
                                    )
                                    morpholine_introduced = True
                                    break
                    else:
                        # Even if not classified as nucleophilic substitution, check for direct evidence
                        # of morpholine introduction
                        if morpholine_reactant and product_has_morpholine:
                            for other_reactant in other_reactants:
                                # Check for any leaving group or halogen
                                has_leaving_group = any(
                                    [
                                        checker.check_fg("Primary halide", other_reactant),
                                        checker.check_fg("Secondary halide", other_reactant),
                                        checker.check_fg("Tertiary halide", other_reactant),
                                        checker.check_fg("Tosylate", other_reactant),
                                        checker.check_fg("Mesylate", other_reactant),
                                        checker.check_fg("Triflate", other_reactant),
                                        "Cl" in other_reactant,
                                        "Br" in other_reactant,
                                        "I" in other_reactant,
                                    ]
                                )

                                if has_leaving_group:
                                    print(
                                        f"Found evidence of morpholine introduction despite reaction type"
                                    )
                                    print(f"Reaction SMILES: {rsmi}")
                                    morpholine_introduced = True
                                    break

                        print(f"Not a nucleophilic substitution reaction at depth {current_depth}")
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return morpholine_introduced
