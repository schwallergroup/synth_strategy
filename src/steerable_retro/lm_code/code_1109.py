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
    This function detects if the synthetic route involves aromatic bromination.
    """
    bromination_detected = False

    def dfs_traverse(node):
        nonlocal bromination_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count bromine atoms in reactants and product
            reactants_br_count = 0
            for r_smiles in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol:
                        for atom in r_mol.GetAtoms():
                            if atom.GetSymbol() == "Br" and atom.GetIsAromatic():
                                reactants_br_count += 1
                except:
                    continue

            # Count bromine atoms in product
            try:
                p_mol = Chem.MolFromSmiles(product_smiles)
                product_br_count = 0
                if p_mol:
                    for atom in p_mol.GetAtoms():
                        if atom.GetSymbol() == "Br":
                            # Check if bromine is attached to aromatic carbon
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetSymbol() == "C" and neighbor.GetIsAromatic():
                                    product_br_count += 1
            except:
                product_br_count = 0

            # If product has more bromine atoms than reactants, it's a bromination
            if product_br_count > reactants_br_count:
                print(f"Aromatic bromination detected in reaction: {rsmi}")
                bromination_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return bromination_detected
