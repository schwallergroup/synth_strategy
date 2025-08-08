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
    This function detects oxazole ring formation in the early stages of synthesis.
    """
    oxazole_formed = False
    formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal oxazole_formed, formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if oxazole is formed in this reaction
            product_mol = Chem.MolFromSmiles(product_smiles)
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

            oxazole_pattern = Chem.MolFromSmarts("[o;r5]1[c;r5][n;r5][c;r5][c;r5]1")

            if product_mol and oxazole_pattern:
                product_has_oxazole = product_mol.HasSubstructMatch(oxazole_pattern)
                reactants_have_oxazole = any(
                    r and r.HasSubstructMatch(oxazole_pattern) for r in reactants_mols if r
                )

                if product_has_oxazole and not reactants_have_oxazole:
                    oxazole_formed = True
                    formation_depth = depth
                    print(f"Oxazole formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Early stage is defined as depth > 3 (closer to starting materials)
    return oxazole_formed and formation_depth > 3
