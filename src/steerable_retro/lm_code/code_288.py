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
    This function detects synthesis routes that involve oxidation of a sulfide to a sulfone.
    """
    has_sulfide_to_sulfone = False

    def dfs_traverse(node, depth=0):
        nonlocal has_sulfide_to_sulfone

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for sulfide to sulfone conversion
            sulfide_pattern = Chem.MolFromSmarts("[#6][S][#6]")
            sulfone_pattern = Chem.MolFromSmarts("[#6][S](=[O])(=[O])[#6]")

            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and sulfone_pattern and product_mol.HasSubstructMatch(sulfone_pattern):
                # Check if sulfide was in reactants
                for r in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and sulfide_pattern and r_mol.HasSubstructMatch(sulfide_pattern):
                        has_sulfide_to_sulfone = True
                        print(f"Found sulfide to sulfone conversion at depth {depth}")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_sulfide_to_sulfone
