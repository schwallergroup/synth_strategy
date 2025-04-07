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
    This function detects a synthetic strategy involving pyrazole formation from a 1,3-diketone and hydrazine.
    """
    pyrazole_formation_found = False

    def dfs_traverse(node):
        nonlocal pyrazole_formation_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains pyrazole
                product_mol = Chem.MolFromSmiles(product)
                pyrazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#7][#6][#6]1")

                if (
                    product_mol
                    and pyrazole_pattern
                    and product_mol.HasSubstructMatch(pyrazole_pattern)
                ):
                    # Check if reactants contain 1,3-diketone and hydrazine
                    diketone_pattern = Chem.MolFromSmarts(
                        "[#6]-[#6](=[O])-[#6]-[#6](=[O])-[#6]"
                    )
                    hydrazine_pattern = Chem.MolFromSmarts("[#7]-[#7]")

                    has_diketone = False
                    has_hydrazine = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if diketone_pattern and reactant_mol.HasSubstructMatch(
                                diketone_pattern
                            ):
                                has_diketone = True
                            if hydrazine_pattern and reactant_mol.HasSubstructMatch(
                                hydrazine_pattern
                            ):
                                has_hydrazine = True

                    if has_diketone and has_hydrazine:
                        print(
                            "Found pyrazole formation from 1,3-diketone and hydrazine"
                        )
                        pyrazole_formation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return pyrazole_formation_found
