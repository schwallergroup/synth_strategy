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
    This function detects if the synthetic route involves halogen exchange transformations.
    """
    halogen_exchange = False
    iodine_pattern = Chem.MolFromSmarts("[#6]-[#53]")
    bromine_pattern = Chem.MolFromSmarts("[#6]-[#35]")
    nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")

    def dfs_traverse(node):
        nonlocal halogen_exchange

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for halogen to nitrile transformation
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol:
                    for r_smi in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smi)
                        if r_mol:
                            # Check if reactant has halogen that's replaced with nitrile in product
                            if (
                                r_mol.HasSubstructMatch(iodine_pattern)
                                or r_mol.HasSubstructMatch(bromine_pattern)
                            ) and product_mol.HasSubstructMatch(nitrile_pattern):
                                print(f"Halogen exchange detected in reaction: {rsmi}")
                                halogen_exchange = True
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return halogen_exchange
