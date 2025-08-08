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
    This function detects a sequence where α-bromination is followed by
    heterocycle formation.
    """
    bromination_reactions = []
    cyclization_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)

                # Check for α-bromination
                if product_mol:
                    alpha_bromo_carbonyl = Chem.MolFromSmarts("[#6][#6](=[#8])[#6][Br]")
                    if product_mol.HasSubstructMatch(alpha_bromo_carbonyl):
                        bromination_reactions.append((depth, node))

                # Check for heterocycle formation (thiazole)
                if product_mol:
                    thiazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#6][#16]1")
                    if product_mol.HasSubstructMatch(thiazole_pattern):
                        cyclization_reactions.append((depth, node))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if bromination is followed by cyclization
    for brom_depth, _ in bromination_reactions:
        for cycl_depth, _ in cyclization_reactions:
            if cycl_depth < brom_depth:  # Lower depth means later in synthesis (closer to product)
                print(
                    f"Detected halogenation-cyclization sequence: bromination at depth {brom_depth}, cyclization at depth {cycl_depth}"
                )
                return True

    return False
