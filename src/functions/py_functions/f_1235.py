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
    This function detects specific functional group interconversions:
    amine to nitro and amide to nitrile.
    """
    amine_to_nitro = False
    amide_to_nitrile = False

    def dfs_traverse(node):
        nonlocal amine_to_nitro, amide_to_nitrile

        if node["type"] == "reaction" and "children" in node:
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if all(reactant_mols) and product_mol:
                # Check for amine to nitro conversion (retrosynthetic direction)
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")

                nitro_in_reactants = any(
                    len(mol.GetSubstructMatches(nitro_pattern)) > 0
                    for mol in reactant_mols
                )
                amine_in_product = (
                    len(product_mol.GetSubstructMatches(amine_pattern)) > 0
                )

                if nitro_in_reactants and amine_in_product:
                    print(f"Amine to nitro conversion detected: {rsmi}")
                    amine_to_nitro = True

                # Check for amide to nitrile conversion (retrosynthetic direction)
                amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")

                nitrile_in_reactants = any(
                    len(mol.GetSubstructMatches(nitrile_pattern)) > 0
                    for mol in reactant_mols
                )
                amide_in_product = (
                    len(product_mol.GetSubstructMatches(amide_pattern)) > 0
                )

                if nitrile_in_reactants and amide_in_product:
                    print(f"Amide to nitrile conversion detected: {rsmi}")
                    amide_to_nitrile = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return amine_to_nitro and amide_to_nitrile
