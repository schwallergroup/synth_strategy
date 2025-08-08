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
    This function detects if the synthetic route involves N-benzyl protection/deprotection.
    """
    benzyl_removed = False

    def dfs_traverse(node):
        nonlocal benzyl_removed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for N-benzyl group in reactants but not in product
            n_benzyl_pattern = Chem.MolFromSmarts("[#7]-[#6]-c1ccccc1")

            reactant_has_nbenzyl = any(
                r is not None and r.HasSubstructMatch(n_benzyl_pattern) for r in reactants
            )
            product_has_nbenzyl = product is not None and product.HasSubstructMatch(
                n_benzyl_pattern
            )

            if reactant_has_nbenzyl and not product_has_nbenzyl:
                print("Detected N-benzyl deprotection")
                benzyl_removed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return benzyl_removed
