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
    This function detects conversion of hydroxyl to chloro groups.
    """
    hydroxyl_to_chloro_found = False

    def dfs_traverse(node):
        nonlocal hydroxyl_to_chloro_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for hydroxyl in reactants
            reactants_mol = Chem.MolFromSmiles(reactants)
            if reactants_mol and reactants_mol.HasSubstructMatch(Chem.MolFromSmarts("c[OH]")):
                # Check for chloro in product at similar position
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c[Cl]")):
                    print("Hydroxyl to chloro conversion detected")
                    hydroxyl_to_chloro_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return hydroxyl_to_chloro_found
