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
    This function detects if the synthetic route involves ester hydrolysis to form a carboxylic acid.
    """
    ester_hydrolysis_found = False

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # First check if this is an ester hydrolysis reaction
            if checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            ):
                print("Detected ester hydrolysis reaction type")

                # Check for ester in reactants
                ester_found = False
                for reactant in reactants:
                    if checker.check_fg("Ester", reactant):
                        print(f"Found ester in reactant: {reactant}")
                        ester_found = True
                        break

                # Check for carboxylic acid in product
                if ester_found and checker.check_fg("Carboxylic acid", product):
                    print(f"Found carboxylic acid in product: {product}")
                    print("Ester hydrolysis confirmed")
                    ester_hydrolysis_found = True

            # Alternative check: look for ester to acid conversion even if reaction type doesn't match
            else:
                print("Not a standard ester hydrolysis reaction, checking functional groups")

                # Check for ester in reactants
                ester_found = False
                for reactant in reactants:
                    if checker.check_fg("Ester", reactant):
                        print(f"Found ester in reactant: {reactant}")
                        ester_found = True
                        break

                # Check for carboxylic acid in product
                if ester_found and checker.check_fg("Carboxylic acid", product):
                    print(f"Found carboxylic acid in product: {product}")
                    print("Ester to acid conversion detected (possible hydrolysis)")
                    ester_hydrolysis_found = True

        # Traverse children
        for child in node.get("children", []):
            if (
                not ester_hydrolysis_found
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: Ester hydrolysis found = {ester_hydrolysis_found}")

    return ester_hydrolysis_found
