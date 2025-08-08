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
    This function detects if the synthesis involves carbamate (Cbz) protection of an amine.
    """
    carbamate_protection_found = False

    # SMARTS patterns
    amine_pattern = Chem.MolFromSmarts("[#7H]")  # NH
    carbamate_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])-[#8]")  # N-C(=O)-O

    def dfs_traverse(node, depth=0):
        nonlocal carbamate_protection_found

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for amine protection with carbamate
            reactant_has_amine = any(
                mol is not None and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
            )
            product_has_carbamate = product_mol is not None and product_mol.HasSubstructMatch(
                carbamate_pattern
            )

            if reactant_has_amine and product_has_carbamate:
                carbamate_protection_found = True
                print(f"Found carbamate protection at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Carbamate protection strategy: {carbamate_protection_found}")
    return carbamate_protection_found
