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
    This function detects if the synthesis route involves sequential functional group
    interconversions, specifically carboxylic acid → acid chloride → amide.
    """
    # Track the sequence of functional group transformations
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Define patterns for functional groups
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                acid_chloride_pattern = Chem.MolFromSmarts("[C](=O)[Cl]")
                amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")

                # Check for transformations
                reactants_have_acid = False
                reactants_have_acid_chloride = False
                product_has_acid_chloride = False
                product_has_amide = False

                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                            reactants_have_acid = True
                        if reactant_mol.HasSubstructMatch(acid_chloride_pattern):
                            reactants_have_acid_chloride = True

                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    if product_mol.HasSubstructMatch(acid_chloride_pattern):
                        product_has_acid_chloride = True
                    if product_mol.HasSubstructMatch(amide_pattern):
                        product_has_amide = True

                # Record transformations with their depth
                if reactants_have_acid and product_has_acid_chloride:
                    transformations.append(("acid_to_acid_chloride", depth))
                    print(f"Carboxylic acid to acid chloride transformation at depth {depth}")

                if reactants_have_acid_chloride and product_has_amide:
                    transformations.append(("acid_chloride_to_amide", depth))
                    print(f"Acid chloride to amide transformation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if we have the sequence acid → acid chloride → amide
    # Sort by depth (higher depth = earlier in synthesis)
    transformations.sort(key=lambda x: x[1], reverse=True)

    # Extract just the transformation types in sequence
    transformation_sequence = [t[0] for t in transformations]

    # Check if our target sequence exists in the correct order
    target_sequence = ["acid_to_acid_chloride", "acid_chloride_to_amide"]

    # Check if target_sequence is a subsequence of transformation_sequence
    result = False
    for i in range(len(transformation_sequence) - len(target_sequence) + 1):
        if transformation_sequence[i : i + len(target_sequence)] == target_sequence:
            result = True
            print(
                "Sequential functional group interconversions detected: acid → acid chloride → amide"
            )
            break

    return result
