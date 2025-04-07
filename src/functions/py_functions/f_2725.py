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
    This function detects late-stage nitrile to carboxylic acid transformation.
    """
    transformation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal transformation_found

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check reactions at depth 0 or 1 (late stage)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile to acid transformation
                nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
                acid_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8;H1,-]")

                product_mol = Chem.MolFromSmiles(product)

                nitrile_found = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitrile_pattern):
                        nitrile_found = True
                        break

                if (
                    nitrile_found
                    and product_mol
                    and product_mol.HasSubstructMatch(acid_pattern)
                ):
                    transformation_found = True
                    print(
                        "Late-stage nitrile to carboxylic acid transformation detected"
                    )

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return transformation_found
