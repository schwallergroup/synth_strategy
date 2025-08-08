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
    This function detects a late-stage SNAr reaction that introduces a piperazine group
    to an aromatic ring, typically replacing a fluorine atom.
    """
    snar_with_piperazine_found = False
    snar_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal snar_with_piperazine_found, snar_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for SNAr with piperazine
            fluoro_aromatic_pattern = Chem.MolFromSmarts("c-F")
            piperazine_pattern = Chem.MolFromSmarts("[N]1CCN([C,c])CC1")
            aromatic_piperazine_pattern = Chem.MolFromSmarts("c-[N]1CCN([C,c])CC1")

            has_fluoro_aromatic = any(
                r.HasSubstructMatch(fluoro_aromatic_pattern) for r in reactants if r
            )
            has_piperazine = any(r.HasSubstructMatch(piperazine_pattern) for r in reactants if r)
            has_aromatic_piperazine = product.HasSubstructMatch(aromatic_piperazine_pattern)

            if has_fluoro_aromatic and has_piperazine and has_aromatic_piperazine:
                snar_with_piperazine_found = True
                snar_depth = depth
                print(f"SNAr with piperazine found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if SNAr with piperazine was found and is late-stage (depth < max_depth/2)
    is_late_stage = snar_depth < max_depth / 2

    if snar_with_piperazine_found and is_late_stage:
        print(f"Found late-stage SNAr with piperazine (depth {snar_depth} out of max {max_depth})")
        return True
    return False
