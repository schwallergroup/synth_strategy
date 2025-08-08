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
    This function detects if an early stage of the synthesis involves
    formation of a sulfonamide bond.
    """
    has_sulfonamide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_sulfonamide_formation

        if node["type"] == "reaction" and depth >= 2:  # Early stage (depth 2 or greater)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if one reactant is a sulfonyl chloride
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])[Cl]")
            amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1]")
            sulfonamide_pattern = Chem.MolFromSmarts("[NX3][SX4](=[OX1])(=[OX1])[#6]")

            has_sulfonyl_chloride = False
            has_amine = False

            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                        has_sulfonyl_chloride = True
                    if mol and mol.HasSubstructMatch(amine_pattern):
                        has_amine = True
                except:
                    continue

            # Check if product has sulfonamide bond
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                has_sulfonamide = product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern)
            except:
                has_sulfonamide = False

            if has_sulfonyl_chloride and has_amine and has_sulfonamide:
                print(f"Detected early-stage sulfonamide formation at depth {depth}")
                has_sulfonamide_formation = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_sulfonamide_formation
