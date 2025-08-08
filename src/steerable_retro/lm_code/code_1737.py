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
    This function detects a synthetic strategy involving ester hydrolysis to form a carboxylic acid.
    """
    found_ester_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ester_hydrolysis

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for ester hydrolysis
            ester_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[O]-[#6]")
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[OH]")

            if product and len(reactants) >= 1:
                # Check if reactant has ester and product has carboxylic acid
                if product.HasSubstructMatch(carboxylic_acid_pattern) and any(
                    r and r.HasSubstructMatch(ester_pattern) for r in reactants if r
                ):
                    print("Found ester hydrolysis")
                    found_ester_hydrolysis = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_ester_hydrolysis
