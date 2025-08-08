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
    Detects ester hydrolysis to carboxylic acid in the synthesis route.
    """
    ester_hydrolysis_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Depth {depth}: Analyzing reaction: {rsmi}")

            # Method 1: Check using reaction type directly
            if checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            ):
                print(f"Depth {depth}: Ester hydrolysis reaction detected directly")
                ester_hydrolysis_detected = True

            # Method 2: Check for ester in reactants and carboxylic acid in product
            else:
                has_ester = False
                for reactant in reactants:
                    if checker.check_fg("Ester", reactant):
                        print(f"Depth {depth}: Ester found in reactant: {reactant}")
                        has_ester = True
                        break

                has_carboxylic_acid = checker.check_fg("Carboxylic acid", product)
                if has_carboxylic_acid:
                    print(f"Depth {depth}: Carboxylic acid found in product: {product}")

                if has_ester and has_carboxylic_acid:
                    # Check for any ester hydrolysis or saponification reaction type
                    if checker.check_reaction(
                        "Ester saponification (methyl deprotection)", rsmi
                    ) or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi):
                        print(f"Depth {depth}: Ester saponification detected")
                        ester_hydrolysis_detected = True
                    else:
                        # If we have an ester in reactants and carboxylic acid in product,
                        # it's likely an ester hydrolysis even if not specifically labeled
                        print(
                            f"Depth {depth}: Potential ester hydrolysis detected based on functional group change"
                        )
                        ester_hydrolysis_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: Ester hydrolysis detected = {ester_hydrolysis_detected}")
    return ester_hydrolysis_detected
