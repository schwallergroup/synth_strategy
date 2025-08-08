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
    Detects a synthetic strategy involving multiple protection steps of a phenol.
    """
    protection_types = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

            # Check for phenol in reactants
            phenol_pattern = Chem.MolFromSmarts("[OH][c]")
            phenol_reactants = any(r and r.HasSubstructMatch(phenol_pattern) for r in reactant_mols)

            # Check for protected phenol patterns in product
            benzyl_pattern = Chem.MolFromSmarts("[c][CH2][O][c]")
            benzoyl_pattern = Chem.MolFromSmarts("[O][C](=[O])[c]")

            if phenol_reactants and product_mol:
                if product_mol.HasSubstructMatch(benzyl_pattern):
                    protection_types.add("benzyl")
                    print(f"Detected benzyl protection at depth {depth}")
                elif product_mol.HasSubstructMatch(benzoyl_pattern):
                    protection_types.add("benzoyl")
                    print(f"Detected benzoyl protection at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have multiple protection types
    multiple_protections = len(protection_types) >= 2

    if multiple_protections:
        print(f"Multiple phenol protection strategy detected: {protection_types}")

    return multiple_protections
