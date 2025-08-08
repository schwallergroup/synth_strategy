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
    This function detects if the final step in the synthesis involves
    urea formation (depth 1 in retrosynthetic analysis).
    """
    late_urea_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_urea_formation

        # Debug the current node
        if node["type"] == "mol":
            print(f"Examining molecule at depth {depth}: {node['smiles'][:30]}...")
        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                print(f"Examining reaction at depth {depth}: {node['metadata']['rsmi'][:30]}...")

        # Check if this is a reaction node at depth 1 (the reaction producing the final product)
        if node["type"] == "reaction" and depth == 1:
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking final reaction: {rsmi[:50]}...")

                # Check if product contains urea or thiourea
                has_urea = checker.check_fg("Urea", product)
                has_thiourea = checker.check_fg("Thiourea", product)

                if has_urea:
                    print("Product contains urea")
                if has_thiourea:
                    print("Product contains thiourea")

                if has_urea or has_thiourea:
                    # Check if the reaction is a known urea formation reaction
                    is_urea_reaction = (
                        checker.check_reaction(
                            "Urea synthesis via isocyanate and primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Urea synthesis via isocyanate and secondary amine", rsmi
                        )
                        or checker.check_reaction("Urea synthesis via isocyanate and diazo", rsmi)
                        or checker.check_reaction(
                            "Urea synthesis via isocyanate and sulfonamide", rsmi
                        )
                        or checker.check_reaction("{urea}", rsmi)
                        or checker.check_reaction("thiourea", rsmi)
                    )

                    # If standard reaction check fails, try alternative detection
                    if not is_urea_reaction:
                        # Check if any reactant contains isocyanate
                        isocyanate_in_reactants = False
                        for reactant in reactants:
                            if checker.check_fg("Isocyanate", reactant):
                                isocyanate_in_reactants = True
                                print("Isocyanate found in reactants")
                                break

                        # If isocyanate is in reactants and urea is in product but not in reactants,
                        # this is likely a urea formation reaction
                        if isocyanate_in_reactants:
                            is_urea_reaction = True
                            print("Detected urea formation via isocyanate reactant")

                    if is_urea_reaction:
                        print("Reaction is a urea formation type")

                        # Verify that urea/thiourea is formed in this reaction (not present in reactants)
                        urea_in_reactants = False
                        for i, reactant in enumerate(reactants):
                            if checker.check_fg("Urea", reactant) or checker.check_fg(
                                "Thiourea", reactant
                            ):
                                print(f"Urea/thiourea found in reactant {i}")
                                urea_in_reactants = True
                                break

                        if not urea_in_reactants:
                            print("Confirmed: Urea/thiourea is formed in this reaction")

                            # Verify this reaction leads to the final product (depth 0)
                            leads_to_final = False
                            for child in node.get("children", []):
                                if child["type"] == "mol" and depth == 1:
                                    leads_to_final = True
                                    break

                            if leads_to_final:
                                print("Detected late-stage urea formation")
                                late_urea_formation = True
                            else:
                                print("Reaction doesn't lead to final product")
                    else:
                        print("Not a known urea formation reaction")
                else:
                    print("Product doesn't contain urea or thiourea")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    if not late_urea_formation:
        print("No late-stage urea formation detected")

    return late_urea_formation
