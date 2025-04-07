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
    This function detects if the synthetic route involves a nitro reduction to amine.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant contains a nitro group
                nitro_pattern = Chem.MolFromSmarts("[#6]c[N+](=[O])[O-]")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(nitro_pattern):
                            # Check if product contains an amine group
                            amine_pattern = Chem.MolFromSmarts("[#6]c[NH2]")
                            prod_mol = Chem.MolFromSmiles(product)
                            if prod_mol and prod_mol.HasSubstructMatch(amine_pattern):
                                print("Nitro reduction to amine detected")
                                nitro_reduction_found = True
                    except:
                        continue

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_found
