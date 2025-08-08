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
    Detects a strategy involving ester hydrolysis as one of the final steps.
    """
    ester_hydrolysis_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester in reactant
            reactant_has_ester = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O][C]")):
                    reactant_has_ester = True
                    break

            # Check for carboxylic acid in product
            if reactant_has_ester:
                products = product.split(".")
                for p in products:
                    prod_mol = Chem.MolFromSmiles(p)
                    if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[OH]")):
                        ester_hydrolysis_depth = depth
                        print(f"Ester hydrolysis detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if ester hydrolysis is in the first two steps (late stage)
    return ester_hydrolysis_depth is not None and ester_hydrolysis_depth <= 1
