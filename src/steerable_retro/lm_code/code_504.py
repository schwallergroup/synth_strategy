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
    This function detects if the synthetic route involves late-stage sulfonamide formation
    as the final step of the synthesis.
    """
    sulfonamide_formed = False

    def dfs_traverse(node):
        nonlocal sulfonamide_formed

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is the final step (depth 0)
            if node.get("depth", 0) == 0:
                # Check for sulfonamide formation
                product_mol = Chem.MolFromSmiles(product)
                sulfonamide_pattern = Chem.MolFromSmarts("[N][S](=O)(=O)[C]")
                if product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern):
                    # Check if sulfonamide wasn't in reactants
                    sulfonamide_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(sulfonamide_pattern):
                            sulfonamide_in_reactants = True
                            break

                    if not sulfonamide_in_reactants:
                        print("Detected late-stage sulfonamide formation")
                        sulfonamide_formed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return sulfonamide_formed
