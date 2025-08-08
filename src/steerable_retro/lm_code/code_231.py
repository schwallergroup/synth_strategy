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
    This function detects late-stage nitro reduction to amine.
    Late stage means it occurs at a low depth (closer to 0).
    """
    nitro_reduction_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro group in reactants
                nitro_in_reactants = False
                for reactant in reactants:
                    if reactant:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
                            if mol.HasSubstructMatch(nitro_pattern):
                                nitro_in_reactants = True

                # Check for amine in product and no nitro
                if product and nitro_in_reactants:
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        amine_pattern = Chem.MolFromSmarts("[NH2]c")
                        nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")

                        if mol.HasSubstructMatch(amine_pattern) and not mol.HasSubstructMatch(
                            nitro_pattern
                        ):
                            nitro_reduction_depth = depth
                            print(f"Nitro reduction detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it occurs at depth 0 or 1
    return nitro_reduction_depth is not None and nitro_reduction_depth <= 1
