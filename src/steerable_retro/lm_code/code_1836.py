#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects if the synthesis involves formation of an oxazole ring.
    """
    # Track if we found oxazole formation
    found_oxazole_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_oxazole_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                products = rsmi.split(">")[-1]

                # Check if product contains oxazole but reactants don't
                product_mol = Chem.MolFromSmiles(products)
                if product_mol:
                    oxazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#6][#8]1")
                    if product_mol.HasSubstructMatch(oxazole_pattern):
                        # Check if reactants don't have oxazole
                        reactants_have_oxazole = False
                        for reactant in reactants.split("."):
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(oxazole_pattern):
                                reactants_have_oxazole = True
                                break

                        if not reactants_have_oxazole:
                            found_oxazole_formation = True
                            print(f"Found oxazole ring formation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_oxazole_formation
