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
    This function detects a synthetic strategy involving amide formation with cyclopropylamine.
    """
    found_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for amide formation with cyclopropylamine
            amide_pattern = Chem.MolFromSmarts("[#6]-[#7]-[#6](=[O])-[#6]")
            cyclopropylamine_pattern = Chem.MolFromSmarts("[NH2]C1CC1")
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[OH]")

            if product and len(reactants) >= 2:
                # Check if product has amide
                if product.HasSubstructMatch(amide_pattern):
                    # Check if reactants include cyclopropylamine and carboxylic acid
                    has_cyclopropylamine = any(
                        r and r.HasSubstructMatch(cyclopropylamine_pattern) for r in reactants if r
                    )
                    has_carboxylic_acid = any(
                        r and r.HasSubstructMatch(carboxylic_acid_pattern) for r in reactants if r
                    )

                    if has_cyclopropylamine and has_carboxylic_acid:
                        print("Found amide formation with cyclopropylamine")
                        found_amide_formation = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_amide_formation
