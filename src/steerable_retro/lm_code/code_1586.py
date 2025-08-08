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
    This function detects if the synthesis includes an ester hydrolysis step.
    """
    has_ester_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal has_ester_hydrolysis

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check if this is an ester hydrolysis reaction using the checker function
            if checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            ):
                print(f"Ester hydrolysis detected at depth {depth}: {rsmi}")
                has_ester_hydrolysis = True

            # Alternative check for ester saponification reactions
            elif checker.check_reaction(
                "Ester saponification (methyl deprotection)", rsmi
            ) or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi):
                print(f"Ester saponification detected at depth {depth}: {rsmi}")
                has_ester_hydrolysis = True

            # If the specific reaction types weren't detected, try a more general approach
            else:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, we're looking for carboxylic acid in reactants and ester in product
                has_carboxylic_acid = any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants if r
                )
                has_alcohol = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Aromatic alcohol", r)
                    for r in reactants
                    if r
                )
                has_ester = checker.check_fg("Ester", product)

                print(
                    f"FG analysis - Carboxylic acid in reactants: {has_carboxylic_acid}, "
                    f"Alcohol in reactants: {has_alcohol}, Ester in product: {has_ester}"
                )

                if has_carboxylic_acid and has_ester:
                    print(
                        f"Potential ester formation (reverse of hydrolysis) detected at depth {depth}: {rsmi}"
                    )
                    if has_alcohol:
                        print(f"Alcohol also present in reactants, confirming ester formation")
                    has_ester_hydrolysis = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_ester_hydrolysis
