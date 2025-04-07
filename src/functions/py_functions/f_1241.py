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
    Detects if the synthesis route uses tert-butyl ester as a protecting group
    for carboxylic acids.
    """
    uses_tert_butyl_protection = False

    def dfs_traverse(node):
        nonlocal uses_tert_butyl_protection

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for tert-butyl ester in reactants
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[O]C(C)(C)C")
                ):
                    uses_tert_butyl_protection = True
                    print(
                        f"Detected tert-butyl ester protection at depth {node['metadata'].get('depth', 'unknown')}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return uses_tert_butyl_protection
