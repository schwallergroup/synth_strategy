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
    Detects if the synthesis route involves multiple ester hydrolysis steps
    as part of a protection/deprotection sequence.
    """
    ester_hydrolysis_count = 0

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester hydrolysis: reactant has ester C(=O)OC pattern and product has carboxylic acid C(=O)O
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    if reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[C](=[O])[O][C]")
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O]")):
                        ester_hydrolysis_count += 1
                        print(
                            f"Detected ester hydrolysis at depth {node['metadata'].get('depth', 'unknown')}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return ester_hydrolysis_count >= 2  # Return True if at least 2 ester hydrolysis steps are found
