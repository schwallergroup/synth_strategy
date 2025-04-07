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
    This function detects O-methylation as a protection strategy for phenols.
    """
    o_methylation_detected = False

    def dfs_traverse(node):
        nonlocal o_methylation_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol in reactants
            phenol_present = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                    if reactant_mol.HasSubstructMatch(phenol_pattern):
                        phenol_present = True
                        break

            # Check for methoxy in product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and phenol_present:
                methoxy_pattern = Chem.MolFromSmarts("[c][O][CH3]")
                if product_mol.HasSubstructMatch(methoxy_pattern):
                    o_methylation_detected = True
                    print("O-methylation detected in reaction:", rsmi)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return o_methylation_detected
