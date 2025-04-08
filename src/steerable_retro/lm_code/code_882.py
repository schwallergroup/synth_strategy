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
    Detects if the synthesis route employs a phenol deprotection strategy,
    where a methoxy group is converted to a hydroxyl group.
    """
    deprotection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal deprotection_detected

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for methoxy in reactants
                methoxy_pattern = Chem.MolFromSmarts("[#6][O][C]")
                # Check for phenol in product
                phenol_pattern = Chem.MolFromSmarts("[#6][O;H1]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if (
                        reactant_mol
                        and product_mol
                        and reactant_mol.HasSubstructMatch(methoxy_pattern)
                        and product_mol.HasSubstructMatch(phenol_pattern)
                        and not reactant_mol.HasSubstructMatch(phenol_pattern)
                    ):
                        print(f"Phenol deprotection detected at depth {depth}")
                        deprotection_detected = True
                        break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return deprotection_detected
