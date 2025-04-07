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
    Detects if the synthesis involves esterification of a carboxylic acid as one of the final steps.
    """
    has_late_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_esterification

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an esterification reaction
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[C](=[O])[O][C]")
            ):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[C](=[O])[OH]")
                    ):
                        # Check if this is a late-stage reaction (depth <= 1)
                        if depth <= 1:
                            has_late_esterification = True
                            print(
                                f"Detected late-stage esterification at depth {depth}"
                            )
                        break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_late_esterification
