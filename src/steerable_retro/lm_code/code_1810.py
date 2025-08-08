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
    This function detects a synthetic strategy involving late-stage epoxide formation
    (in the final or penultimate step of the synthesis).
    """
    epoxide_formation_found = False
    epoxide_formation_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal epoxide_formation_found, epoxide_formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for epoxide formation
            product_mol = Chem.MolFromSmiles(product_smiles)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

            # Epoxide pattern
            epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")

            if product_mol.HasSubstructMatch(epoxide_pattern):
                # Check if epoxide is newly formed (not present in reactants)
                epoxide_in_reactants = any(
                    r.HasSubstructMatch(epoxide_pattern) for r in reactant_mols if r is not None
                )

                if not epoxide_in_reactants:
                    epoxide_formation_found = True
                    epoxide_formation_depth = min(epoxide_formation_depth, depth)
                    print(f"Epoxide formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if epoxide formation is late-stage (depth 0 or 1)
    late_stage = epoxide_formation_found and epoxide_formation_depth <= 1

    if late_stage:
        print("Late-stage epoxide formation strategy detected")

    return late_stage
