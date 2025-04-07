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
    This function detects nitro reduction to amine in the synthesis route.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro group in reactants
                nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
                amine_pattern = Chem.MolFromSmarts("[#7;H2]")

                # Check if any reactant has nitro group
                reactant_has_nitro = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if (
                        reactant_mol
                        and nitro_pattern
                        and reactant_mol.HasSubstructMatch(nitro_pattern)
                    ):
                        reactant_has_nitro = True
                        break

                # Check if product has new amine group
                product_mol = Chem.MolFromSmiles(product)
                if (
                    reactant_has_nitro
                    and product_mol
                    and amine_pattern
                    and product_mol.HasSubstructMatch(amine_pattern)
                ):
                    nitro_reduction_found = True
                    print("Nitro reduction detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_found
