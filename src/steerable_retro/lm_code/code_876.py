#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    Detects if the synthesis involves late-stage sulfonamide formation as the final step.

    This function checks if the final reaction step (depth 1) is a sulfonamide formation
    reaction where a sulfonyl chloride reacts with an amine to form a sulfonamide.
    """
    found_sulfonamide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_sulfonamide_formation

        # Skip if we've already found what we're looking for
        if found_sulfonamide_formation:
            return

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            print(f"Checking reaction at depth {depth}: {node['metadata'].get('rsmi', 'No RSMI')}")

            # The final reaction step in retrosynthesis is at depth 1
            if depth == 1:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Checking final reaction step (depth 1): {rsmi}")

                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if this is a sulfonamide synthesis reaction using predefined reactions
                    is_sulfonamide_reaction = False

                    # Check for known sulfonamide formation reactions
                    sulfonamide_reactions = [
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                        "Schotten-Baumann to ester",  # This might include some sulfonamide formations
                    ]

                    for reaction_name in sulfonamide_reactions:
                        if checker.check_reaction(reaction_name, rsmi):
                            print(f"Detected sulfonamide formation via {reaction_name}")
                            is_sulfonamide_reaction = True
                            break

                    # Fallback: Check reactants and products manually
                    if not is_sulfonamide_reaction:
                        print("Checking reactants and product manually")
                        has_sulfonyl_chloride = any(
                            checker.check_fg("Sulfonyl halide", r) for r in reactants
                        )

                        # Check for various amine types
                        has_amine = any(
                            checker.check_fg(amine_type, r)
                            for r in reactants
                            for amine_type in [
                                "Primary amine",
                                "Secondary amine",
                                "Aniline",
                                "Tertiary amine",
                            ]
                        )

                        # Check for sulfonamide in product
                        has_sulfonamide = checker.check_fg("Sulfonamide", product)

                        print(f"Has sulfonyl halide: {has_sulfonyl_chloride}")
                        print(f"Has amine: {has_amine}")
                        print(f"Has sulfonamide in product: {has_sulfonamide}")

                        if has_sulfonyl_chloride and has_amine and has_sulfonamide:
                            print("Detected sulfonamide formation from reactants and product")
                            is_sulfonamide_reaction = True

                    if is_sulfonamide_reaction:
                        print("Found late-stage sulfonamide formation")
                        found_sulfonamide_formation = True
                except Exception as e:
                    print(f"Error processing reaction node: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            # In retrosynthesis, depth increases as we move from products to reactants
            next_depth = depth + 1
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)
    return found_sulfonamide_formation
