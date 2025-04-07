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
    transformations from carboxylic acid to ester to alcohol to mesylate.
    """
    # Track if we've seen each functional group transformation
    acid_to_ester = False
    ester_to_alcohol = False
    alcohol_to_mesylate = False

    # SMARTS patterns for functional groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("[#6][C](=[O])[O;H]")
    methyl_ester_pattern = Chem.MolFromSmarts("[#6][C](=[O])[O][C;H3]")
    alcohol_pattern = Chem.MolFromSmarts("[#6][C;H2][O;H]")
    mesylate_pattern = Chem.MolFromSmarts("[#6][C;H2][O][S](=[O])(=[O])[C;H3]")

    def dfs_traverse(node):
        nonlocal acid_to_ester, ester_to_alcohol, alcohol_to_mesylate

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactants_mol and product_mol:
                # Check for acid to ester transformation
                if (
                    reactants_mol.HasSubstructMatch(carboxylic_acid_pattern)
                    and product_mol.HasSubstructMatch(methyl_ester_pattern)
                    and not reactants_mol.HasSubstructMatch(methyl_ester_pattern)
                ):
                    acid_to_ester = True
                    print("Carboxylic acid to ester transformation detected")

                # Check for ester to alcohol transformation
                if (
                    reactants_mol.HasSubstructMatch(methyl_ester_pattern)
                    and product_mol.HasSubstructMatch(alcohol_pattern)
                    and not reactants_mol.HasSubstructMatch(alcohol_pattern)
                ):
                    ester_to_alcohol = True
                    print("Ester to alcohol transformation detected")

                # Check for alcohol to mesylate transformation
                if (
                    reactants_mol.HasSubstructMatch(alcohol_pattern)
                    and product_mol.HasSubstructMatch(mesylate_pattern)
                    and not reactants_mol.HasSubstructMatch(mesylate_pattern)
                ):
                    alcohol_to_mesylate = True
                    print("Alcohol to mesylate transformation detected")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the complete sequence is present
    sequence_present = acid_to_ester and ester_to_alcohol and alcohol_to_mesylate

    if sequence_present:
        print("Complete carboxylic acid → ester → alcohol → mesylate sequence detected")

    return sequence_present
