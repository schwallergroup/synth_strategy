#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects the use of trifluoroacetamide as a protecting group for amines.
    It looks for the formation and subsequent removal of trifluoroacetamide groups.
    """
    trifluoroacetamide_pattern = Chem.MolFromSmarts("[NH]C(=O)C(F)(F)F")
    trifluoroacetamide_formations = 0

    def dfs_traverse(node):
        nonlocal trifluoroacetamide_formations

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for trifluoroacetamide formation
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(trifluoroacetamide_pattern):
                # Check if trifluoroacetamide was not in reactants
                trifluoroacetamide_in_reactants = False
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol and reactant_mol.HasSubstructMatch(trifluoroacetamide_pattern):
                        trifluoroacetamide_in_reactants = True
                        break

                if not trifluoroacetamide_in_reactants:
                    print("Trifluoroacetamide formation detected")
                    trifluoroacetamide_formations += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return trifluoroacetamide_formations >= 1
