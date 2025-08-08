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
    Detects a synthetic strategy involving pyrazole heterocycle disconnection.
    """
    has_pyrazole_disconnection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_pyrazole_disconnection

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(r for r in reactants):
                pyrazole_pattern = Chem.MolFromSmarts("[#7]1-[#7]=[#6]-[#6]-[#6]1")

                if product.HasSubstructMatch(pyrazole_pattern) and any(
                    r.HasSubstructMatch(pyrazole_pattern) for r in reactants
                ):
                    print(f"Found pyrazole disconnection at depth {depth}")
                    has_pyrazole_disconnection = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if has_pyrazole_disconnection:
        print("Detected pyrazole-containing synthesis strategy")
    else:
        print("Pyrazole-containing synthesis strategy not detected")

    return has_pyrazole_disconnection
