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
    This function detects a synthetic strategy involving addition of a
    trifluoromethyl group to a carbonyl.
    """
    trifluoromethyl_addition_found = False

    def dfs_traverse(node):
        nonlocal trifluoromethyl_addition_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ketone in reactants
            ketone_pattern = Chem.MolFromSmarts("[C](=[O])[#6]")

            for r_smi in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r_smi)
                if r_mol and r_mol.HasSubstructMatch(ketone_pattern):
                    # Check for CF3-containing alcohol in product
                    cf3_alcohol_pattern = Chem.MolFromSmarts("[C]([OH])([#6])[C]([F])([F])[F]")
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(cf3_alcohol_pattern):
                        print("Found trifluoromethyl addition to carbonyl")
                        trifluoromethyl_addition_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return trifluoromethyl_addition_found
