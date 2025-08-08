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
    Detects if the synthesis route includes a debromination step (C-Br to C-H conversion).
    """
    found_debromination = False

    def dfs_traverse(node):
        nonlocal found_debromination

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for C-Br bond in reactants
            c_br_pattern = Chem.MolFromSmarts("[#6][Br]")

            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(c_br_pattern):
                    # Check if the bromine is absent in the product
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        # If product has fewer bromines than reactant, it's likely a debromination
                        reactant_br_count = len(reactant_mol.GetSubstructMatches(c_br_pattern))
                        product_br_count = len(product_mol.GetSubstructMatches(c_br_pattern))

                        if product_br_count < reactant_br_count:
                            found_debromination = True
                            print("Found debromination step")
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_debromination
