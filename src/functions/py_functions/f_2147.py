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
    This function detects the strategy of forming an oxadiazole heterocycle.
    """
    oxadiazole_formed = False

    def dfs_traverse(node):
        nonlocal oxadiazole_formed

        if node["type"] == "reaction":
            # Check if this is a reaction node with metadata
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride or similar in reactants
                acid_chloride_pattern = Chem.MolFromSmarts("[#6](=[O])[#17]")
                has_acid_chloride = False

                for reactant in reactants:
                    if reactant:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(acid_chloride_pattern):
                                has_acid_chloride = True
                                break
                        except:
                            continue

                # Check for oxadiazole in product
                if has_acid_chloride:
                    oxadiazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#7][#6][#8]1")
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            oxadiazole_pattern
                        ):
                            print("Detected oxadiazole formation")
                            oxadiazole_formed = True
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return oxadiazole_formed
