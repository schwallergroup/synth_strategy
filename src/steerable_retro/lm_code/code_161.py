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


def main(route):
    """
    This function detects if the synthesis involves carbamate protection of an amine.
    """
    carbamate_pattern = Chem.MolFromSmarts("[#7][#6](=[#8])[#8]")
    amine_to_carbamate = False

    def dfs_traverse(node, depth=0):
        nonlocal amine_to_carbamate

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amine in reactants
                amine_in_reactants = any("N" in r and not "[N+]" in r for r in reactants)

                # Check for carbamate in product
                product_mol = Chem.MolFromSmiles(product)
                carbamate_in_product = product_mol and product_mol.HasSubstructMatch(
                    carbamate_pattern
                )

                if amine_in_reactants and carbamate_in_product:
                    amine_to_carbamate = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Carbamate protection strategy: {amine_to_carbamate}")
    return amine_to_carbamate
