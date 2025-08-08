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
    Detects if the synthesis involves amide formation from an ester.
    """
    amide_formation_detected = False

    def dfs_traverse(node):
        nonlocal amide_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester pattern in reactants
            ester_pattern = Chem.MolFromSmarts("[#6][#8]C(=O)[#6]")

            # Check for amide pattern in product
            amide_pattern = Chem.MolFromSmarts("[#7]C(=O)[#6]")

            # Check if any reactant has ester group
            ester_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(ester_pattern):
                    ester_in_reactants = True
                    break

            # Check if product has amide group
            product_mol = Chem.MolFromSmiles(product)
            amide_in_product = product_mol and product_mol.HasSubstructMatch(amide_pattern)

            if ester_in_reactants and amide_in_product:
                print("Amide formation from ester detected")
                amide_formation_detected = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_detected
