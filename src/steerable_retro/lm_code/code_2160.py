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
    Detects a linear synthesis strategy with sequential functionalization:
    phenol → ether → aldehyde → alcohol → alkyne → coupled product
    """
    # Track functional group transformations
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to molecules
            reactant_mols = [
                Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
            ]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and reactant_mols:
                # Define patterns for functional groups
                phenol_pattern = Chem.MolFromSmarts("[OX2H][cX3]:[c]")
                ether_pattern = Chem.MolFromSmarts("[OX2]([cX3])[CX4]")
                aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
                alcohol_pattern = Chem.MolFromSmarts("[OX2H][CX4]")
                terminal_alkyne_pattern = Chem.MolFromSmarts("[CX2]#[CH]")
                internal_alkyne_pattern = Chem.MolFromSmarts("[CX2]#[CX2]")

                # Check for transformations
                if any(
                    r.HasSubstructMatch(phenol_pattern) for r in reactant_mols
                ) and product_mol.HasSubstructMatch(ether_pattern):
                    transformations.append(("phenol_to_ether", depth))
                    print(f"Found phenol to ether transformation at depth {depth}")

                if any(
                    r.HasSubstructMatch(ether_pattern) for r in reactant_mols
                ) and product_mol.HasSubstructMatch(aldehyde_pattern):
                    transformations.append(("ether_to_aldehyde", depth))
                    print(f"Found ether to aldehyde transformation at depth {depth}")

                if any(
                    r.HasSubstructMatch(aldehyde_pattern) for r in reactant_mols
                ) and product_mol.HasSubstructMatch(alcohol_pattern):
                    transformations.append(("aldehyde_to_alcohol", depth))
                    print(f"Found aldehyde to alcohol transformation at depth {depth}")

                if any(
                    r.HasSubstructMatch(alcohol_pattern) for r in reactant_mols
                ) and product_mol.HasSubstructMatch(terminal_alkyne_pattern):
                    transformations.append(("alcohol_to_terminal_alkyne", depth))
                    print(f"Found alcohol to terminal alkyne transformation at depth {depth}")

                if any(
                    r.HasSubstructMatch(terminal_alkyne_pattern) for r in reactant_mols
                ) and product_mol.HasSubstructMatch(internal_alkyne_pattern):
                    transformations.append(("terminal_to_internal_alkyne", depth))
                    print(f"Found terminal to internal alkyne transformation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 3 sequential transformations
    if len(transformations) >= 3:
        # Sort by depth to check sequence
        transformations.sort(
            key=lambda x: x[1], reverse=True
        )  # Higher depth = earlier in synthesis
        print(f"Found sequential transformations: {[t[0] for t in transformations]}")
        return True
    return False
