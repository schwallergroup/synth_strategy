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
    This function detects if the synthetic route involves sequential functionalization
    of an aromatic ring (at least 3 sequential transformations on the same ring).
    """
    # Track transformations on aromatic rings
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check if reaction modifies an aromatic ring
                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Patterns for common aromatic transformations
                    patterns = [
                        ("nitration", Chem.MolFromSmarts("c[N+](=[O])[O-]")),
                        ("reduction", Chem.MolFromSmarts("c[NH2]")),
                        ("alkylation", Chem.MolFromSmarts("c[O][C]")),
                        ("halogenation", Chem.MolFromSmarts("c[F,Cl,Br,I]")),
                        ("amination", Chem.MolFromSmarts("c[N]")),
                    ]

                    for name, pattern in patterns:
                        if product_mol.HasSubstructMatch(pattern):
                            transformations.append((depth, name))
                            print(f"Found {name} at depth {depth}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Sort transformations by depth (ascending order - early to late in synthesis)
    transformations.sort(key=lambda x: x[0], reverse=True)

    # Check if we have at least 3 sequential transformations
    return len(transformations) >= 3
