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
    This function detects if the synthetic route involves protection and/or deprotection of a phenol group.
    Specifically looking for MOM (methoxymethyl) protection.
    """
    phenol_pattern = Chem.MolFromSmarts("Oc1ccccc1")
    protected_phenol_pattern = Chem.MolFromSmarts("COCOc1ccccc1")
    has_phenol_protection = False

    def dfs_traverse(node):
        nonlocal has_phenol_protection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check for protection (phenol -> protected phenol)
                    if (
                        reactants_mol
                        and product_mol
                        and reactants_mol.HasSubstructMatch(phenol_pattern)
                        and product_mol.HasSubstructMatch(protected_phenol_pattern)
                    ):
                        print(f"Detected phenol protection in reaction: {rsmi}")
                        has_phenol_protection = True

                    # Check for deprotection (protected phenol -> phenol)
                    if (
                        reactants_mol
                        and product_mol
                        and reactants_mol.HasSubstructMatch(protected_phenol_pattern)
                        and product_mol.HasSubstructMatch(phenol_pattern)
                    ):
                        print(f"Detected phenol deprotection in reaction: {rsmi}")
                        has_phenol_protection = True
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_phenol_protection
