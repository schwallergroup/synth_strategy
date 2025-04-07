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
    This function detects if a synthetic route contains benzylic functionalization
    (specifically tolyl to benzyl bromide conversion) as a key step.
    """
    found_benzylic_functionalization = False

    def dfs_traverse(node):
        nonlocal found_benzylic_functionalization

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for tolyl to benzyl bromide conversion
            tolyl_pattern = Chem.MolFromSmarts("[c][C]")
            benzyl_bromide_pattern = Chem.MolFromSmarts("[c][C][Br]")

            has_tolyl = False
            for reactant in reactants:
                if reactant.strip():
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(tolyl_pattern):
                            has_tolyl = True
                            break
                    except:
                        continue

            has_benzyl_bromide = False
            if product.strip():
                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol and mol.HasSubstructMatch(benzyl_bromide_pattern):
                        has_benzyl_bromide = True
                except:
                    pass

            # If transformation is tolyl to benzyl bromide
            if has_tolyl and has_benzyl_bromide:
                found_benzylic_functionalization = True
                print(f"Found benzylic functionalization: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Benzylic functionalization strategy detected: {found_benzylic_functionalization}")
    return found_benzylic_functionalization
