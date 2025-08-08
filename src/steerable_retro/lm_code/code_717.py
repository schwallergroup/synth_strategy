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
    This function detects the use of geminal dinitrile (malononitrile) as a key building block
    and its manipulation throughout the synthesis.
    """
    # Initialize flags
    has_geminal_dinitrile = False
    manipulates_geminal_dinitrile = False

    def dfs_traverse(node):
        nonlocal has_geminal_dinitrile, manipulates_geminal_dinitrile

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for geminal dinitrile group
                geminal_dinitrile_pattern = Chem.MolFromSmarts("[#6]([C]#[N])([C]#[N])")
                if mol.HasSubstructMatch(geminal_dinitrile_pattern):
                    has_geminal_dinitrile = True
                    print(f"Found geminal dinitrile group in molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
            product = Chem.MolFromSmiles(product_str) if product_str else None

            if product and reactants:
                geminal_dinitrile_pattern = Chem.MolFromSmarts("[#6]([C]#[N])([C]#[N])")

                # Check if reaction involves manipulation of geminal dinitrile
                reactant_has_geminal = any(
                    mol and mol.HasSubstructMatch(geminal_dinitrile_pattern) for mol in reactants
                )
                product_has_geminal = product and product.HasSubstructMatch(
                    geminal_dinitrile_pattern
                )

                if reactant_has_geminal and product_has_geminal:
                    manipulates_geminal_dinitrile = True
                    print(f"Found reaction manipulating geminal dinitrile: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Traverse the route
    dfs_traverse(route)

    # Check if geminal dinitrile is present and manipulated
    strategy_detected = has_geminal_dinitrile and manipulates_geminal_dinitrile

    if strategy_detected:
        print("Detected geminal dinitrile strategy with manipulation throughout synthesis")
    else:
        print("Did not detect complete geminal dinitrile strategy")
        print(
            f"Has geminal dinitrile: {has_geminal_dinitrile}, Manipulates geminal dinitrile: {manipulates_geminal_dinitrile}"
        )

    return strategy_detected
