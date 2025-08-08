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
    Detects if the synthesis includes an amide to nitrile conversion step.
    """
    has_conversion = False

    def dfs_traverse(node, depth=0):
        nonlocal has_conversion

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Extract reactants and products
            try:
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]

                reactant_smiles_list = reactants_part.split(".")
                product_smiles = products_part

                # Check if any reactant has an amide group and the product has a nitrile group
                for reactant_smiles in reactant_smiles_list:
                    has_amide = (
                        checker.check_fg("Primary amide", reactant_smiles)
                        or checker.check_fg("Secondary amide", reactant_smiles)
                        or checker.check_fg("Tertiary amide", reactant_smiles)
                    )

                    has_nitrile = checker.check_fg("Nitrile", product_smiles)

                    if has_amide:
                        print(f"Found amide in reactant: {reactant_smiles}")
                    if has_nitrile:
                        print(f"Found nitrile in product: {product_smiles}")

                    if has_amide and has_nitrile:
                        # Check for known reaction types that could convert amide to nitrile
                        known_reactions = [
                            "Dehydration",
                            "Oxidation of amide to carboxylic acid",
                            "Reduction of amides to amines",  # Could be part of a multi-step process
                            "Hydrogenolysis of amides/imides/carbamates",  # Could be relevant
                        ]

                        for reaction_type in known_reactions:
                            if checker.check_reaction(reaction_type, rsmi):
                                print(f"Found matching reaction type: {reaction_type}")
                                has_conversion = True
                                return

                        # If no specific reaction type is found but we have amide → nitrile transformation
                        # This is a fallback check for when the specific reaction type isn't in our list
                        print(
                            f"Found amide to nitrile conversion without specific reaction type match"
                        )
                        has_conversion = True
                        return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {has_conversion}")
    return has_conversion
