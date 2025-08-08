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
    This function detects amide formation from acid chloride and amine.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acid chloride in reactants
            acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")

            # Check for amide in product
            amide_pattern = Chem.MolFromSmarts("C(=O)N")

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("N")

            # Check if reactants have acid chloride and amine, and product has amide
            reactants_have_acid_chloride = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(acid_chloride_pattern)
                for r in reactants
                if r
            )
            reactants_have_amine = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern) for r in reactants if r
            )
            product_has_amide = (
                Chem.MolFromSmiles(product).HasSubstructMatch(amide_pattern) if product else False
            )

            if reactants_have_acid_chloride and reactants_have_amine and product_has_amide:
                amide_formation_found = True
                print("Amide formation from acid chloride detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Amide formation from acid chloride detected: {amide_formation_found}")
    return amide_formation_found
