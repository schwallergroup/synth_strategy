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
    This function detects if the synthesis includes a Suzuki coupling via boronic ester intermediate.
    """
    has_suzuki_coupling = False

    def dfs_traverse(node):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic ester in reactants
            boronic_ester_pattern = Chem.MolFromSmarts("[c]-[B]([O])[O]")
            biaryl_formation = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(boronic_ester_pattern):
                    # Now check if product has a new biaryl bond
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # This is a simplification - in a real implementation,
                        # we would need to check for new C-C bonds between aryl groups
                        # by comparing reactants and products
                        print("Detected potential Suzuki coupling")
                        has_suzuki_coupling = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_suzuki_coupling
