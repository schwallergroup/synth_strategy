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
    Detects the formation of a quinazoline ring system through condensation
    of an aldehyde and an amine.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde in reactants
                aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[NX3;H2]")
                # Check for quinazoline in product
                quinazoline_pattern = Chem.MolFromSmarts("c1ncnc2ccccc12")

                has_aldehyde = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(aldehyde_pattern):
                            has_aldehyde = True
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(quinazoline_pattern):
                    if has_aldehyde and has_amine:
                        found_pattern = True
                        print(f"Found quinazoline formation via condensation at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_pattern
