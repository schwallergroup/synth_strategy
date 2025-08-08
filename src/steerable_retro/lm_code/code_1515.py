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
    This function detects a sequential functional group transformation pattern:
    Br → CN → COOH → amide
    """
    # Track transformations
    transformations = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]

            if product_mol and reactant_mols:
                # Check for specific transformations

                # Br → CN transformation
                if any(
                    mol.HasSubstructMatch(Chem.MolFromSmarts("c[Br]")) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("cC#N")):
                    transformations.append("Br_to_CN")
                    print("Detected Br → CN transformation")

                # CN → COOH transformation
                if any(
                    mol.HasSubstructMatch(Chem.MolFromSmarts("cC#N")) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("cC(=O)O")):
                    transformations.append("CN_to_COOH")
                    print("Detected CN → COOH transformation")

                # COOH → amide transformation
                if any(
                    mol.HasSubstructMatch(Chem.MolFromSmarts("cC(=O)O")) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("cC(=O)N")):
                    transformations.append("COOH_to_amide")
                    print("Detected COOH → amide transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if we have the complete sequence
    has_sequence = all(t in transformations for t in ["Br_to_CN", "CN_to_COOH", "COOH_to_amide"])

    print(f"Sequential functional group transformation strategy detected: {has_sequence}")
    print(f"Transformations found: {transformations}")

    return has_sequence
