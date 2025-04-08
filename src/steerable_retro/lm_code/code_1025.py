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
    This function detects a three-step sequence involving:
    1. Introduction of sulfur via nucleophilic substitution
    2. Oxidation to sulfoxide
    3. Further oxidation to sulfone
    """
    # Initialize tracking variables
    reactions_sequence = []

    # SMARTS patterns
    thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")
    sulfoxide_pattern = Chem.MolFromSmarts("[#6]-[#16](=[#8])-[#6]")
    sulfone_pattern = Chem.MolFromSmarts("[#6]-[#16](=[#8])(=[#8])-[#6]")

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Determine reaction type
                reaction_type = None

                # Check for C-S bond formation (thioether formation)
                if product_mol.HasSubstructMatch(thioether_pattern) and not any(
                    mol.HasSubstructMatch(thioether_pattern) for mol in reactant_mols
                ):
                    reaction_type = "thioether_formation"

                # Check for thioether to sulfoxide oxidation
                elif (
                    any(mol.HasSubstructMatch(thioether_pattern) for mol in reactant_mols)
                    and product_mol.HasSubstructMatch(sulfoxide_pattern)
                    and not any(mol.HasSubstructMatch(sulfoxide_pattern) for mol in reactant_mols)
                ):
                    reaction_type = "sulfoxide_formation"

                # Check for sulfoxide to sulfone oxidation
                elif (
                    any(mol.HasSubstructMatch(sulfoxide_pattern) for mol in reactant_mols)
                    and product_mol.HasSubstructMatch(sulfone_pattern)
                    and not any(mol.HasSubstructMatch(sulfone_pattern) for mol in reactant_mols)
                ):
                    reaction_type = "sulfone_formation"

                if reaction_type:
                    reactions_sequence.append((depth, reaction_type))
                    print(f"Detected {reaction_type} at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (to get chronological order in synthesis direction)
    reactions_sequence.sort(key=lambda x: x[0], reverse=True)

    # Extract just the reaction types in sequence
    reaction_types = [r[1] for r in reactions_sequence]

    # Check if we have the complete sequence in the correct order
    target_sequence = ["thioether_formation", "sulfoxide_formation", "sulfone_formation"]

    # Check if target_sequence is a subsequence of reaction_types
    is_subsequence = False
    if len(reaction_types) >= len(target_sequence):
        for i in range(len(reaction_types) - len(target_sequence) + 1):
            if reaction_types[i : i + len(target_sequence)] == target_sequence:
                is_subsequence = True
                break

    print(f"Three-step sulfur oxidation sequence detected: {is_subsequence}")
    print(f"Reaction sequence: {reaction_types}")

    return is_subsequence
