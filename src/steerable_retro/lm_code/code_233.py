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
    This function detects a late-stage amide to nitrile conversion strategy.
    Late stage is defined as occurring at depth 0 or 1.
    """
    found_conversion = False

    def dfs_traverse(node):
        nonlocal found_conversion

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            depth = node.get("depth", -1)

            # Check for amide to nitrile conversion
            amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
            nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")

            amide_found = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(amide_pattern):
                    amide_found = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            nitrile_found = product_mol and product_mol.HasSubstructMatch(nitrile_pattern)

            # Check if this is a late-stage conversion (depth 0 or 1)
            if amide_found and nitrile_found and depth <= 1:
                found_conversion = True
                print(f"Found late-stage amide to nitrile conversion at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage amide to nitrile conversion strategy detected: {found_conversion}")
    return found_conversion
