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
    This function detects if the synthetic route involves benzimidazole ring formation.
    """
    benzimidazole_pattern = Chem.MolFromSmarts("c1nc2ccccc2[nH]1")
    benzimidazole_formed = False

    def dfs_traverse(node):
        nonlocal benzimidazole_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants don't have benzimidazole but product does
                reactants_have_benzimidazole = False
                for r_smiles in reactants_smiles:
                    try:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol and r_mol.HasSubstructMatch(benzimidazole_pattern):
                            reactants_have_benzimidazole = True
                            break
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    product_has_benzimidazole = (
                        product_mol
                        and product_mol.HasSubstructMatch(benzimidazole_pattern)
                    )

                    if not reactants_have_benzimidazole and product_has_benzimidazole:
                        print("Benzimidazole formation detected")
                        benzimidazole_formed = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return benzimidazole_formed
