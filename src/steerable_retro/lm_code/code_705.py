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
    This function detects if thiourea formation occurs in the late stage of synthesis
    (low depth in the synthetic tree).
    """
    found_late_thiourea = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_thiourea

        if node["type"] == "reaction" and depth <= 1:  # Consider depth 0 or 1 as late stage
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for thiourea formation
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                thiourea_pattern = Chem.MolFromSmarts("NC(=S)N")
                if product_mol.HasSubstructMatch(thiourea_pattern):
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")
                            if reactant_mol.HasSubstructMatch(isothiocyanate_pattern):
                                found_late_thiourea = True
                                print(f"Found thiourea formation at depth {depth}")
                                break

        # Recursively traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_late_thiourea
