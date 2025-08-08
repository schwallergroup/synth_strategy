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
    This function detects Wittig-type olefination of an aldehyde.
    """
    has_wittig_olefination = False

    def dfs_traverse(node, depth=0):
        nonlocal has_wittig_olefination

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                products = rsmi.split(">")[-1]

                # Check for Wittig olefination
                if "." in reactants and "P+" in reactants:  # Phosphonium salt present
                    reactant_mol = Chem.MolFromSmiles(reactants)
                    product_mol = Chem.MolFromSmiles(products)

                    if reactant_mol and product_mol:
                        # Aldehyde pattern
                        aldehyde_pattern = Chem.MolFromSmarts("[#6]=O")
                        # Alkene pattern
                        alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")

                        if "P+" in reactants and "[CH" in reactants and "O=" in reactants:
                            if product_mol.HasSubstructMatch(alkene_pattern):
                                print(f"Detected Wittig-type olefination at depth {depth}")
                                has_wittig_olefination = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_wittig_olefination
