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
    This function detects a synthetic strategy involving nitrile hydrolysis to
    carboxylic acid or derivative.
    """
    nitrile_hydrolysis_detected = False

    def dfs_traverse(node):
        nonlocal nitrile_hydrolysis_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitrile in reactants
            nitrile_in_reactants = False
            for r_smi in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smi)
                    if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("C#N")):
                        nitrile_in_reactants = True
                        break
                except:
                    continue

            # Check for carboxylic acid in product
            if nitrile_in_reactants:
                try:
                    p_mol = Chem.MolFromSmiles(product_smiles)
                    if p_mol and p_mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O")):
                        nitrile_hydrolysis_detected = True
                        print("Nitrile hydrolysis detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nitrile_hydrolysis_detected
