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
    This function detects if the synthetic route involves tert-butyl ester deprotection.
    """
    tert_butyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC(C)(C)C")
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    deprotection_detected = False

    def dfs_traverse(node):
        nonlocal deprotection_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if any reactant contains tert-butyl ester
            reactant_has_tert_butyl = False
            for r_smiles in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol and r_mol.HasSubstructMatch(tert_butyl_ester_pattern):
                        reactant_has_tert_butyl = True
                        break
                except:
                    continue

            # Check if product contains carboxylic acid
            try:
                p_mol = Chem.MolFromSmiles(product_smiles)
                product_has_carboxylic_acid = p_mol and p_mol.HasSubstructMatch(
                    carboxylic_acid_pattern
                )
            except:
                product_has_carboxylic_acid = False

            # If reactant has tert-butyl ester and product has carboxylic acid, it's a deprotection
            if reactant_has_tert_butyl and product_has_carboxylic_acid:
                print(f"tert-butyl ester deprotection detected in reaction: {rsmi}")
                deprotection_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return deprotection_detected
