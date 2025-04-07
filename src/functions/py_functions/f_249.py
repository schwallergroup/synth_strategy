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
    This function detects sequential ester hydrolysis steps in the synthesis.
    """
    ester_hydrolysis_count = 0

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester hydrolysis pattern
            try:
                reactant_mol = Chem.MolFromSmiles(
                    reactants[0]
                )  # Assuming first reactant
                product_mol = Chem.MolFromSmiles(product)

                # SMARTS patterns for ester and carboxylic acid
                ester_pattern = Chem.MolFromSmarts("C(=O)OC")
                acid_pattern = Chem.MolFromSmarts("C(=O)O")

                # Check if reactant has ester and product has acid
                if (
                    reactant_mol
                    and reactant_mol.HasSubstructMatch(ester_pattern)
                    and product_mol
                    and product_mol.HasSubstructMatch(acid_pattern)
                ):
                    print(
                        f"Detected ester hydrolysis at depth {node['metadata'].get('depth', -1)}"
                    )
                    ester_hydrolysis_count += 1
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return (
        ester_hydrolysis_count >= 2
    )  # Return True if at least 2 ester hydrolysis steps
