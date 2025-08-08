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
    This function detects if the synthesis involves halogen displacement reactions.
    """
    halogen_displacement_count = 0

    def dfs_traverse(node):
        nonlocal halogen_displacement_count

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for halogen in reactants but not in product (or fewer halogens in product)
            reactant_mol = Chem.MolFromSmiles(reactants_part)
            product_mol = Chem.MolFromSmiles(product_part)

            if reactant_mol and product_mol:
                # Count halogens in reactants and products
                halogen_pattern = Chem.MolFromSmarts("[#9,#17,#35,#53]")  # F, Cl, Br, I
                reactant_halogens = len(reactant_mol.GetSubstructMatches(halogen_pattern))
                product_halogens = len(product_mol.GetSubstructMatches(halogen_pattern))

                if reactant_halogens > product_halogens:
                    halogen_displacement_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    result = halogen_displacement_count > 0

    print(f"Halogen displacement reactions: {result} (count: {halogen_displacement_count})")
    return result
