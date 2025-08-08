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
    Detects O-methylation protection of hydroxyl groups in the synthesis route.
    """
    protection_found = False

    def dfs_traverse(node):
        nonlocal protection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for OH to OCH3 transformation
            reactant_mol = None
            for r in reactants:
                if "[OH:" in r:
                    reactant_mol = Chem.MolFromSmiles(r)
                    break

            if reactant_mol and "[O:" in product and "[CH3:" in product:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check if this is an O-methylation
                    for atom in reactant_mol.GetAtoms():
                        if atom.GetSymbol() == "O" and atom.GetTotalNumHs() > 0:
                            protection_found = True
                            print("O-methylation protection detected")
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return protection_found
