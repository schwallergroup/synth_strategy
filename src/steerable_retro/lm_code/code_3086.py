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
    Detects if the synthesis route involves formation of a sulfonamide functional group.
    """
    has_sulfonamide_formation = False

    def dfs_traverse(node):
        nonlocal has_sulfonamide_formation

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for sulfonamide formation
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#17]")
            sulfonamide_pattern = Chem.MolFromSmarts("[#7][#16](=[#8])(=[#8])")

            # Check if any reactant contains sulfonyl chloride
            has_sulfonyl_chloride = any(
                Chem.MolFromSmiles(reactant) is not None
                and Chem.MolFromSmiles(reactant).HasSubstructMatch(sulfonyl_chloride_pattern)
                for reactant in reactants_smiles
                if reactant
            )

            # Check if product contains sulfonamide
            product_mol = Chem.MolFromSmiles(product_smiles)
            has_sulfonamide = product_mol is not None and product_mol.HasSubstructMatch(
                sulfonamide_pattern
            )

            if has_sulfonyl_chloride and has_sulfonamide:
                has_sulfonamide_formation = True
                print(f"Detected sulfonamide formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Sulfonamide-containing strategy detected: {has_sulfonamide_formation}")
    return has_sulfonamide_formation
