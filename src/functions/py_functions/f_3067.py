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
    Detects a linear synthesis with sequential functional group transformations:
    Cl → CN → COOH → CH2OH followed by etherification
    """
    # Track if we've found each transformation
    found_transformations = {
        "chloride_to_nitrile": False,
        "nitrile_to_acid": False,
        "acid_to_alcohol": False,
        "alcohol_to_ether": False,
    }

    # Track the depth at which each transformation occurs
    transformation_depths = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for chloride to nitrile transformation
            if any(
                mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4][Cl]"))
                for mol in reactants
            ) and product.HasSubstructMatch(Chem.MolFromSmarts("[CX4][CX2]#[NX1]")):
                found_transformations["chloride_to_nitrile"] = True
                transformation_depths["chloride_to_nitrile"] = depth
                print(f"Found chloride to nitrile transformation at depth {depth}")

            # Check for nitrile to carboxylic acid transformation
            if any(
                mol.HasSubstructMatch(Chem.MolFromSmarts("[CX2]#[NX1]"))
                for mol in reactants
            ) and product.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")):
                found_transformations["nitrile_to_acid"] = True
                transformation_depths["nitrile_to_acid"] = depth
                print(
                    f"Found nitrile to carboxylic acid transformation at depth {depth}"
                )

            # Check for carboxylic acid to alcohol transformation
            if any(
                mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]"))
                for mol in reactants
            ) and product.HasSubstructMatch(Chem.MolFromSmarts("[CX4][OX2H]")):
                found_transformations["acid_to_alcohol"] = True
                transformation_depths["acid_to_alcohol"] = depth
                print(
                    f"Found carboxylic acid to alcohol transformation at depth {depth}"
                )

            # Check for alcohol to ether transformation
            if any(
                mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4][OX2H]"))
                for mol in reactants
            ) and product.HasSubstructMatch(Chem.MolFromSmarts("[CX4][OX2][c]")):
                found_transformations["alcohol_to_ether"] = True
                transformation_depths["alcohol_to_ether"] = depth
                print(f"Found alcohol to ether transformation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found all transformations in the correct sequence
    if all(found_transformations.values()):
        # Verify the sequence (higher depth = earlier in synthesis)
        depths = transformation_depths
        if (
            depths["chloride_to_nitrile"]
            > depths["nitrile_to_acid"]
            > depths["acid_to_alcohol"]
            > depths["alcohol_to_ether"]
        ):
            print("Found complete sequential functional group transformation strategy")
            return True

    return False
