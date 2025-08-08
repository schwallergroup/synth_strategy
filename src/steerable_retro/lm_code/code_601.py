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
    This function detects if the synthesis involves a late-stage reductive amination
    (conversion of a ketone to a secondary amine in the final steps).
    """
    # Initialize tracking variables
    found_reductive_amination = False
    depth_of_amination = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_reductive_amination, depth_of_amination

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for reductive amination using the checker function
                if checker.check_reaction("Reductive amination with ketone", rsmi):
                    print(f"Found reductive amination with ketone at depth {depth}")
                    found_reductive_amination = True
                    depth_of_amination = min(depth_of_amination, depth)
                elif checker.check_reaction("Reductive amination with aldehyde", rsmi):
                    print(f"Found reductive amination with aldehyde at depth {depth}")
                    found_reductive_amination = True
                    depth_of_amination = min(depth_of_amination, depth)

                # Fallback check if the specific reaction types aren't detected
                if not found_reductive_amination:
                    product_mol = Chem.MolFromSmiles(product)
                    has_ketone = False
                    has_amine = False
                    has_secondary_amine_in_product = False

                    if product_mol and checker.check_fg("Secondary amine", product):
                        has_secondary_amine_in_product = True
                        print(f"Product has secondary amine at depth {depth}")

                    for reactant in reactants:
                        if checker.check_fg("Ketone", reactant):
                            has_ketone = True
                            print(f"Reactant has ketone at depth {depth}")
                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Aniline", reactant
                        ):
                            has_amine = True
                            print(f"Reactant has amine at depth {depth}")

                    if has_ketone and has_amine and has_secondary_amine_in_product:
                        print(f"Detected likely reductive amination at depth {depth}")
                        found_reductive_amination = True
                        depth_of_amination = min(depth_of_amination, depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it's in the first few steps of the synthesis (low depth in retrosynthetic direction)
    max_depth = 2  # Allow reactions up to depth 2 to be considered late-stage
    result = found_reductive_amination and depth_of_amination <= max_depth

    print(f"Found reductive amination: {found_reductive_amination}, at depth: {depth_of_amination}")
    print(f"Is late-stage reductive amination: {result}")

    return result
