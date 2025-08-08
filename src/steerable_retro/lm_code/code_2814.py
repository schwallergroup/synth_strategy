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
    Detects if the synthesis involves conversion of an ester to an alcohol.
    """
    ester_to_alcohol_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_to_alcohol_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                # Check for ester in reactants
                ester_pattern = Chem.MolFromSmarts("[#6][O][C]=[O]")
                has_ester = any(
                    mol.HasSubstructMatch(ester_pattern) for mol in reactant_mols if mol
                )

                # Check for alcohol in product
                alcohol_pattern = Chem.MolFromSmarts("[#6][OH]")
                has_alcohol = product_mol and product_mol.HasSubstructMatch(alcohol_pattern)

                if has_ester and has_alcohol:
                    print(f"Detected ester to alcohol conversion at depth {depth}")
                    ester_to_alcohol_detected = True
            except:
                print("Error processing reaction SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return ester_to_alcohol_detected
