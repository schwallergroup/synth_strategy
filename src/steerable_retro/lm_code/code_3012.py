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
    This function detects a synthetic strategy involving triazole ring formation
    in the middle stage of synthesis.
    """
    triazole_formed = False

    def dfs_traverse(node):
        nonlocal triazole_formed

        if node["type"] == "reaction":
            # Check if this is a middle-stage reaction (depth around 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert SMILES to molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product:
                    # Check if product contains a triazole
                    triazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#7][#6][#6]1")
                    if product.HasSubstructMatch(triazole_pattern):
                        # Check if reactants don't have triazole
                        reactants_have_triazole = False
                        for r in reactants:
                            if r and r.HasSubstructMatch(triazole_pattern):
                                reactants_have_triazole = True
                                break

                        if not reactants_have_triazole:
                            print("Detected triazole formation")
                            triazole_formed = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return triazole_formed
