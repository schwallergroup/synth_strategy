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
    This function detects if the synthetic route maintains a lactam core structure
    throughout the synthesis.
    """
    lactam_reactions_count = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal lactam_reactions_count, total_reactions

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                total_reactions += 1
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for lactam pattern in product
                lactam_pattern = Chem.MolFromSmarts("[#6]-[#7]-[#6](=[#8])-[#6]")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(lactam_pattern):
                        print(f"Found lactam core in reaction product: {product}")
                        lactam_reactions_count += 1
                except:
                    print(f"Error processing SMILES: {product}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if lactam core is present in at least 80% of reactions
    return total_reactions > 0 and (lactam_reactions_count / total_reactions) >= 0.8
