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
    This function detects a synthetic strategy involving thiazole ring formation.
    """
    thiazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal thiazole_formation_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for thiazole in product
            product_mol = Chem.MolFromSmiles(product)
            thiazole_pattern = Chem.MolFromSmarts("c1scnc1")

            # Check if thiazole exists in product but not in reactants
            if product_mol and product_mol.HasSubstructMatch(thiazole_pattern):
                thiazole_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(thiazole_pattern):
                        thiazole_in_reactants = True
                        break

                if not thiazole_in_reactants:
                    thiazole_formation_detected = True
                    print(f"Detected thiazole formation at depth {node.get('depth', 'unknown')}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Thiazole formation detected: {thiazole_formation_detected}")
    return thiazole_formation_detected
