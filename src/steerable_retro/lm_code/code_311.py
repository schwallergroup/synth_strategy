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
    This function detects a late-stage epoxide opening reaction that forms a secondary alcohol.
    Late stage is defined as occurring at depth 0 or 1 in the synthesis tree.
    """
    epoxide_opening_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal epoxide_opening_detected

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for epoxide in reactants
                epoxide_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")
                # Check for alcohol in product
                alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8H]")

                epoxide_in_reactants = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(epoxide_pattern):
                            epoxide_in_reactants = True
                            break
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    alcohol_in_product = product_mol and product_mol.HasSubstructMatch(
                        alcohol_pattern
                    )
                except:
                    alcohol_in_product = False

                if epoxide_in_reactants and alcohol_in_product:
                    print(f"Detected epoxide opening at depth {depth}")
                    epoxide_opening_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return epoxide_opening_detected
