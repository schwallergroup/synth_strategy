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
    Detects if the synthesis includes an amine to nitrile functional group transformation.
    """
    amine_to_nitrile = False

    def dfs_traverse(node, depth=0):
        nonlocal amine_to_nitrile

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                reaction_smiles = node["metadata"]["rsmi"]
                reactants_part = reaction_smiles.split(">")[0]
                products_part = reaction_smiles.split(">")[-1]

                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[N;H2]")
                reactant_mol = Chem.MolFromSmiles(reactants_part)

                # Check for nitrile in products
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                product_mol = Chem.MolFromSmiles(products_part)

                if (
                    reactant_mol
                    and product_mol
                    and reactant_mol.HasSubstructMatch(amine_pattern)
                    and product_mol.HasSubstructMatch(nitrile_pattern)
                ):
                    amine_to_nitrile = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Amine to nitrile transformation: {amine_to_nitrile}")

    return amine_to_nitrile
