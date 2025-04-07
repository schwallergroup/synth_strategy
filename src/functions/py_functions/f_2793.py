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
    Detects a specific sequence of functional group transformations:
    tertiary alcohol → alkene → aldehyde → alkene → alcohol/ketone
    """
    # Track functional group transformations
    transformations = []

    # SMARTS patterns for functional groups
    tertiary_alcohol_pattern = Chem.MolFromSmarts("[#6]-[#6](-[OH])-[#6]")
    alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    aldehyde_pattern = Chem.MolFromSmarts("[#6]-[#6]=[O]")

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [
                    Chem.MolFromSmiles(r)
                    for r in reactants_smiles
                    if Chem.MolFromSmiles(r)
                ]

                # Identify functional groups in reactants and products
                reactant_groups = []
                for mol in reactant_mols:
                    if mol:
                        if mol.HasSubstructMatch(tertiary_alcohol_pattern):
                            reactant_groups.append("tertiary_alcohol")
                        if mol.HasSubstructMatch(alkene_pattern):
                            reactant_groups.append("alkene")
                        if mol.HasSubstructMatch(aldehyde_pattern):
                            reactant_groups.append("aldehyde")

                product_groups = []
                if product_mol:
                    if product_mol.HasSubstructMatch(tertiary_alcohol_pattern):
                        product_groups.append("tertiary_alcohol")
                    if product_mol.HasSubstructMatch(alkene_pattern):
                        product_groups.append("alkene")
                    if product_mol.HasSubstructMatch(aldehyde_pattern):
                        product_groups.append("aldehyde")

                # Record transformation
                if reactant_groups and product_groups:
                    transformation = (
                        f"{'+'.join(reactant_groups)} → {'+'.join(product_groups)}"
                    )
                    transformations.append((depth, transformation))

            except:
                print("Error processing SMILES in reaction at depth", depth)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth
    transformations.sort(key=lambda x: x[0])

    # Check for the specific sequence
    sequence_present = False

    # Print all transformations for debugging
    print("Functional group transformations (depth, transformation):")
    for depth, transformation in transformations:
        print(f"  {depth}: {transformation}")

    # Look for patterns that match our sequence
    alcohol_to_alkene = False
    alkene_to_aldehyde = False
    aldehyde_to_alkene = False

    for i in range(len(transformations) - 1):
        current = transformations[i][1]
        next_trans = transformations[i + 1][1]

        if (
            "tertiary_alcohol" in current
            and "alkene" in next_trans
            and "tertiary_alcohol" not in next_trans
        ):
            alcohol_to_alkene = True

        if (
            "alkene" in current
            and "aldehyde" in next_trans
            and "alkene" not in next_trans
        ):
            alkene_to_aldehyde = True

        if (
            "aldehyde" in current
            and "alkene" in next_trans
            and "aldehyde" not in next_trans
        ):
            aldehyde_to_alkene = True

    sequence_present = alcohol_to_alkene and alkene_to_aldehyde and aldehyde_to_alkene

    print(f"Functional group sequence detection results:")
    print(f"  Alcohol to alkene: {alcohol_to_alkene}")
    print(f"  Alkene to aldehyde: {alkene_to_aldehyde}")
    print(f"  Aldehyde to alkene: {aldehyde_to_alkene}")
    print(f"  Complete sequence present: {sequence_present}")

    return sequence_present
