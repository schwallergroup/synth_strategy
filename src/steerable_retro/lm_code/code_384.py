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
    This function detects a linear synthesis strategy that includes
    an ester hydrolysis step in preparation for amide formation.
    """
    has_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal has_ester_hydrolysis

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ester hydrolysis pattern
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            # Look for ester in reactants
            ester_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")

            # Look for carboxylic acid in product
            acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")

            has_ester = any(mol and mol.HasSubstructMatch(ester_pattern) for mol in reactant_mols)
            has_acid = product_mol and product_mol.HasSubstructMatch(acid_pattern)

            if has_ester and has_acid:
                print("Detected ester hydrolysis step")
                has_ester_hydrolysis = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_ester_hydrolysis
