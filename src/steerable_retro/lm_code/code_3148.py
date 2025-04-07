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
    Detects the use of amide bond formation as a key strategy in the synthesis.
    """
    amide_formations = 0

    def dfs_traverse(node):
        nonlocal amide_formations

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains an amide bond
                amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    # Count amide bonds in product
                    product_amide_count = len(product_mol.GetSubstructMatches(amide_pattern))

                    # Count amide bonds in reactants
                    reactant_amide_count = 0
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                reactant_amide_count += len(mol.GetSubstructMatches(amide_pattern))
                        except:
                            continue

                    # If product has more amide bonds than reactants combined, amide formation occurred
                    if product_amide_count > reactant_amide_count:
                        print(
                            f"Detected amide bond formation: {reactant_amide_count} → {product_amide_count}"
                        )
                        amide_formations += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total amide bond formations: {amide_formations}")
    return amide_formations >= 1
