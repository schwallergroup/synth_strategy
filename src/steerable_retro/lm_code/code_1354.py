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
    This function detects a synthetic strategy involving late-stage introduction
    of a morpholine-containing fragment.
    """
    morpholine_introduction = False

    def dfs_traverse(node, depth=0):
        nonlocal morpholine_introduction

        if node["type"] == "reaction" and depth <= 1:  # Late-stage reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for morpholine in reactants
                morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1")

                has_morpholine_reactant = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(morpholine_pattern):
                            has_morpholine_reactant = True
                            break
                    except:
                        continue

                # Check if product also has morpholine
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        has_morpholine_reactant
                        and product_mol
                        and product_mol.HasSubstructMatch(morpholine_pattern)
                    ):
                        print(f"Morpholine introduction detected at depth {depth}")
                        morpholine_introduction = True
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return morpholine_introduction
