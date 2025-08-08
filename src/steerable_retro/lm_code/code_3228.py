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
    This function detects if the route uses a linear fragment assembly strategy
    with at least 3 fragments combined sequentially.
    """
    fragment_coupling_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if this is a coupling reaction (2+ reactants)
            if len(reactants) >= 2:
                # Look for common coupling patterns
                product = rsmi.split(">")[-1]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Check for amide coupling
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                    # Check for C-N bond formation
                    cn_bond_pattern = Chem.MolFromSmarts("[C][N]")

                    if product_mol.HasSubstructMatch(
                        amide_pattern
                    ) or product_mol.HasSubstructMatch(cn_bond_pattern):
                        fragment_coupling_depths.append(depth)
                        print(f"Found fragment coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have at least 2 coupling reactions (to combine 3+ fragments)
    # and that they occur in a sequential manner (different depths)
    if len(fragment_coupling_depths) >= 2 and len(set(fragment_coupling_depths)) >= 2:
        print(f"Found linear fragment assembly with {len(fragment_coupling_depths)} couplings")
        return True
    return False
