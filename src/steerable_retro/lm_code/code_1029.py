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
    This function detects if a protection step (specifically Boc protection)
    occurs in the middle of the synthesis.
    """
    boc_protection_depth = None
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for Boc protection
                boc_pattern = Chem.MolFromSmarts("[CH3]C([CH3])([CH3])[O][C](=[O])[N]")

                if product.HasSubstructMatch(boc_pattern) and not any(
                    r.HasSubstructMatch(boc_pattern) for r in reactants
                ):
                    boc_protection_depth = depth
                    print(f"Detected Boc protection at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if Boc protection occurs in the middle of the synthesis
    # In retrosynthetic direction, middle means not at depth 0 and not at max_depth
    is_mid_synthesis = (
        boc_protection_depth is not None
        and boc_protection_depth > 0
        and boc_protection_depth < max_depth
    )

    print(f"Mid-synthesis protection detected: {is_mid_synthesis}")
    return is_mid_synthesis
