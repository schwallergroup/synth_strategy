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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects if the synthetic route involves protection of carboxylic acid as ester.
    """
    carboxylic_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal carboxylic_protection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Method 1: Check if this is a protection reaction directly
                if checker.check_reaction("Protection of carboxylic acid", rsmi):
                    print(f"Carboxylic acid protection detected via reaction type check: {rsmi}")
                    carboxylic_protection = True
                    return

                # Method 2: Check for esterification reactions that could be protection
                if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                    print(f"Potential carboxylic acid protection via esterification: {rsmi}")
                    carboxylic_protection = True
                    return

                # Method 3: Check functional group transformation manually
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and checker.check_fg("Ester", product_smiles):
                    for r_smi in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", r_smi):
                            print(
                                f"Carboxylic acid protection detected via FG transformation: {rsmi}"
                            )
                            print(f"  Reactant with carboxylic acid: {r_smi}")
                            print(f"  Product with ester: {product_smiles}")
                            carboxylic_protection = True
                            return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    print("Starting traversal to find carboxylic acid protection...")
    dfs_traverse(route)
    print(f"Carboxylic acid protection found: {carboxylic_protection}")
    return carboxylic_protection
