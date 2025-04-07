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
    This function detects if the synthetic route involves a sequence of
    nitro reduction to amine followed by halogenation.
    """
    # Track if we've seen each transformation
    nitro_to_amine = False
    amine_to_halogen = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine, amine_to_halogen

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            try:
                reactant_mol = Chem.MolFromSmiles(reactants_part)
                product_mol = Chem.MolFromSmiles(products_part)

                if reactant_mol and product_mol:
                    # Check for nitro to amine transformation
                    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

                    if reactant_mol.HasSubstructMatch(
                        nitro_pattern
                    ) and product_mol.HasSubstructMatch(amine_pattern):
                        print("Nitro to amine transformation detected")
                        nitro_to_amine = True

                    # Check for amine to halogen transformation
                    halogen_pattern = Chem.MolFromSmarts("[#6]-[#9,#17,#35,#53]")

                    if reactant_mol.HasSubstructMatch(
                        amine_pattern
                    ) and product_mol.HasSubstructMatch(halogen_pattern):
                        print("Amine to halogen transformation detected")
                        amine_to_halogen = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_to_amine and amine_to_halogen
