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
    This function detects if the synthesis involves formation of aryl ethers.
    """
    aryl_ether_formations = 0

    def dfs_traverse(node):
        nonlocal aryl_ether_formations

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for aryl ether formation
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[F,Cl,Br,I]")
            phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")
            aryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")

            reactants_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactants_mol and product_mol:
                if (
                    reactants_mol.HasSubstructMatch(aryl_halide_pattern)
                    or reactants_mol.HasSubstructMatch(phenol_pattern)
                ) and product_mol.HasSubstructMatch(aryl_ether_pattern):
                    aryl_ether_formations += 1
                    print(f"Aryl ether formation found, total: {aryl_ether_formations}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if at least one aryl ether formation was found
    if aryl_ether_formations >= 1:
        print(f"Aryl ether formation strategy detected with {aryl_ether_formations} formations")
        return True
    return False
