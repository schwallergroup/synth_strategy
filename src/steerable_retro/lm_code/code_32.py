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
    This function detects if aromatic halogens (Cl, F) are preserved throughout the synthesis.
    """
    halogen_preserved = True
    reactions_checked = 0

    def dfs_traverse(node):
        nonlocal halogen_preserved, reactions_checked

        if node["type"] == "reaction":
            reactions_checked += 1

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            product_mol = Chem.MolFromSmiles(product_smiles)

            if not product_mol:
                print("Warning: Could not parse product molecule")
                return

            # Check for aromatic chlorine and fluorine
            cl_pattern = Chem.MolFromSmarts("c[Cl]")
            f_pattern = Chem.MolFromSmarts("c[F]")

            if not (
                product_mol.HasSubstructMatch(cl_pattern)
                and product_mol.HasSubstructMatch(f_pattern)
            ):
                print("Aromatic halogen not preserved in a reaction")
                halogen_preserved = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Only consider valid if we checked at least one reaction
    if reactions_checked == 0:
        return False

    print(f"Aromatic halogen preservation: {halogen_preserved}")
    return halogen_preserved
