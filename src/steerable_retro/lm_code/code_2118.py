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
    This function detects if the synthetic route involves SNAr chemistry for C-N bond formation.
    """
    snar_detected = False

    def dfs_traverse(node):
        nonlocal snar_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aromatic C-F bond in reactants
            aromatic_cf_pattern = Chem.MolFromSmarts("c-F")

            # Check for aromatic C-N bond in product where F was
            aromatic_cn_pattern = Chem.MolFromSmarts("c-[#7]")

            # Check if any reactant has an aromatic C-F bond
            for reactant_smiles in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                if reactant_mol and reactant_mol.HasSubstructMatch(aromatic_cf_pattern):
                    # Check if product has a new aromatic C-N bond
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(aromatic_cn_pattern):
                        # This is a simplified check - ideally we'd confirm it's the same carbon
                        snar_detected = True
                        print(f"Detected potential SNAr reaction: {rsmi}")
                        break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"SNAr for C-N bond formation detected: {snar_detected}")
    return snar_detected
