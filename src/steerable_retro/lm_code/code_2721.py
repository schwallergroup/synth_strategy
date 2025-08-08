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
    This function detects if the synthetic route employs a late-stage sulfonamide formation
    as the final step (depth 0) or penultimate step (depth 1).
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # Extract depth with proper error handling
            depth_match = re.search(
                r"Depth: (\d+)", node.get("metadata", {}).get("ID", "Depth: -1")
            )
            depth = int(depth_match.group(1)) if depth_match else -1

            if depth <= 1:  # Final or penultimate step (late-stage)
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if product contains sulfonamide
                product_has_sulfonamide = checker.check_fg("Sulfonamide", product_smiles)

                if product_has_sulfonamide:
                    print(f"Product contains sulfonamide at depth {depth}")

                    # Check if sulfonamide is formed in this step (not present in reactants)
                    sulfonamide_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Sulfonamide", reactant):
                            sulfonamide_in_reactants = True
                            print(f"Sulfonamide already present in reactant: {reactant}")
                            break

                    # Check for sulfonyl chloride in reactants
                    sulfonyl_chloride_present = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Sulfonyl halide", reactant):
                            sulfonyl_chloride_present = True
                            print(f"Found sulfonyl chloride in reactant: {reactant}")
                            break

                    # Check for amine in reactants
                    amine_present = False
                    for reactant in reactants_smiles:
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        ):
                            amine_present = True
                            print(f"Found amine in reactant: {reactant}")
                            break

                    # Check if this is a sulfonamide formation reaction
                    is_sulfonamide_reaction = checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    ) or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    )

                    if is_sulfonamide_reaction:
                        print(f"Reaction is a sulfonamide formation reaction")

                    if not sulfonamide_in_reactants and (
                        is_sulfonamide_reaction or (sulfonyl_chloride_present and amine_present)
                    ):
                        print(f"Detected late-stage sulfonamide formation at depth {depth}")
                        result = True
                else:
                    print(f"Product does not contain sulfonamide at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: {result}")
    return result
