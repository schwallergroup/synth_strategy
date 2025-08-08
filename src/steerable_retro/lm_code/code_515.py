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
    This function detects a synthetic strategy where a nitro group is introduced
    to an aromatic ring in an intermediate step and then reduced to an amine
    in the final step of the synthesis.
    """
    # Track if we found nitration and reduction steps
    nitration_found = False
    nitro_reduction_found = False

    # Track the depth of these reactions
    nitration_depth = -1
    reduction_depth = -1

    print("Starting analysis of synthesis route")

    def dfs_traverse(node, depth=0):
        nonlocal nitration_found, nitro_reduction_found, nitration_depth, reduction_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and products
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for nitration reaction
                is_nitration = False
                nitration_reactions = [
                    "Aromatic nitration with HNO3",
                    "Aromatic nitration with NO3 salt",
                    "Aromatic nitration with NO2 salt",
                    "Aromatic nitration with alkyl NO2",
                ]

                for nitration_type in nitration_reactions:
                    if checker.check_reaction(nitration_type, rsmi):
                        is_nitration = True
                        print(f"Found nitration reaction: {nitration_type}")
                        break

                # Check if nitro group is being added (not present in reactants but present in product)
                reactants_have_nitro = any(
                    checker.check_fg("Nitro group", r) for r in reactants_smiles
                )
                product_has_nitro = checker.check_fg("Nitro group", product_smiles)

                if is_nitration or (not reactants_have_nitro and product_has_nitro):
                    nitration_found = True
                    nitration_depth = depth
                    print(f"Confirmed nitration reaction at depth {depth}")

                # Check for nitro reduction (nitro to amine)
                is_nitro_reduction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )

                # If not a standard reduction, check for functional group changes
                if not is_nitro_reduction:
                    reactants_have_nitro = any(
                        checker.check_fg("Nitro group", r) for r in reactants_smiles
                    )
                    product_has_nitro = checker.check_fg("Nitro group", product_smiles)
                    product_has_amine = checker.check_fg(
                        "Primary amine", product_smiles
                    ) or checker.check_fg("Aniline", product_smiles)

                    if reactants_have_nitro and not product_has_nitro and product_has_amine:
                        is_nitro_reduction = True
                        print("Detected nitro reduction based on functional group changes")

                if is_nitro_reduction:
                    nitro_reduction_found = True
                    reduction_depth = depth
                    print(f"Confirmed nitro reduction at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Analysis complete: nitration_found={nitration_found} at depth {nitration_depth}, "
        f"nitro_reduction_found={nitro_reduction_found} at depth {reduction_depth}"
    )

    # Check if we found both nitration and reduction in the correct order
    # Remember: lower depth is later in the synthesis (closer to final product)
    if nitration_found and nitro_reduction_found and reduction_depth < nitration_depth:
        print("Late-stage nitro reduction strategy detected")
        return True

    return False
