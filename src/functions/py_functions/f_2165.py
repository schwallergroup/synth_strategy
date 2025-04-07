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
    Detects if the synthesis route involves connecting heterocycles via an oxygen linker.
    """
    found_oxygen_linker = False

    def dfs_traverse(node):
        nonlocal found_oxygen_linker

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for formation of c-O-C pattern connecting heterocycles
            oxygen_linker_pattern = Chem.MolFromSmarts("c-[#8]-[#6]")
            benzothiazole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nc(s2)")

            product_mol = Chem.MolFromSmiles(product)
            if not product_mol:
                return

            if product_mol.HasSubstructMatch(
                oxygen_linker_pattern
            ) and product_mol.HasSubstructMatch(benzothiazole_pattern):
                # Check if this pattern was formed in this reaction
                pattern_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue
                    if reactant_mol.HasSubstructMatch(
                        oxygen_linker_pattern
                    ) and reactant_mol.HasSubstructMatch(benzothiazole_pattern):
                        pattern_in_reactants = True
                        break

                if not pattern_in_reactants:
                    found_oxygen_linker = True
                    print(
                        f"Found oxygen linker formation at depth: {node.get('depth', 'unknown')}"
                    )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return found_oxygen_linker
