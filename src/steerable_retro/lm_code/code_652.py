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
    This function detects the conversion of a hydroxyl group to a chloro group.
    """
    oh_to_cl_conversion_detected = False

    def dfs_traverse(node):
        nonlocal oh_to_cl_conversion_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for OH in reactants and Cl in the same position in product
            for reactant_smiles in reactants_smiles:
                if "[OH]" in reactant_smiles or "OH" in reactant_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if reactant_mol and product_mol:
                        # This is a simplified check - in a real implementation,
                        # you would need to map atoms between reactants and products
                        if "[Cl]" in product_smiles or "Cl" in product_smiles:
                            print(f"OH to Cl conversion detected in reaction: {rsmi}")
                            oh_to_cl_conversion_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return oh_to_cl_conversion_detected
