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


def main(route, min_count=2):
    """
    Detects if the synthetic route involves multiple amide bond formations.
    """
    # Count amide formations
    amide_formation_count = 0

    def is_amide_formation(reaction_smiles):
        """Check if a reaction forms an amide bond"""
        # Split into reactants and product
        parts = reaction_smiles.split(">")
        if len(parts) < 3:
            return False

        reactants = parts[0].split(".")
        product = parts[2]

        # Check for carboxylic acid in reactants
        carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8;H1]")

        # Check for amine in reactants
        amine_pattern = Chem.MolFromSmarts("[#7;!$(NC=O)]")

        # Check for amide in product
        amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")

        has_acid = False
        has_amine = False

        for reactant in reactants:
            try:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(carboxylic_acid_pattern):
                        has_acid = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True
            except:
                continue

        has_amide = False
        try:
            prod_mol = Chem.MolFromSmiles(product)
            if prod_mol and prod_mol.HasSubstructMatch(amide_pattern):
                has_amide = True
        except:
            pass

        return (has_acid or has_amine) and has_amide

    def dfs_traverse(node):
        nonlocal amide_formation_count

        # Check if this is a reaction node
        if node.get("type") == "reaction":
            # Get reaction SMILES from metadata
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check if this is an amide formation
                if is_amide_formation(rsmi):
                    amide_formation_count += 1
                    print(f"Found amide formation (count: {amide_formation_count})")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total amide formations: {amide_formation_count}")
    return amide_formation_count >= min_count
