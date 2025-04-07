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
    This function detects a strategy involving sequential oxidation of sulfur:
    thioether → sulfoxide → sulfone, following introduction of sulfur via mesylate activation.
    """
    # Initialize tracking variables
    has_alcohol_to_mesylate = False
    has_mesylate_to_thioether = False
    has_thioether_to_sulfoxide = False
    has_sulfoxide_to_sulfone = False

    # SMARTS patterns for functional groups
    alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
    mesylate_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#16](=[#8])(=[#8])-[#6]")
    thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")
    sulfoxide_pattern = Chem.MolFromSmarts("[#6]-[#16](=[#8])-[#6]")
    sulfone_pattern = Chem.MolFromSmarts("[#6]-[#16](=[#8])(=[#8])-[#6]")

    def dfs_traverse(node):
        nonlocal has_alcohol_to_mesylate, has_mesylate_to_thioether
        nonlocal has_thioether_to_sulfoxide, has_sulfoxide_to_sulfone

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check for alcohol to mesylate transformation
                if any(
                    mol.HasSubstructMatch(alcohol_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(mesylate_pattern):
                    print("Detected alcohol to mesylate transformation")
                    has_alcohol_to_mesylate = True

                # Check for mesylate to thioether transformation
                if (
                    any(
                        mol.HasSubstructMatch(mesylate_pattern) for mol in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(thioether_pattern)
                    and not product_mol.HasSubstructMatch(mesylate_pattern)
                ):
                    print("Detected mesylate to thioether transformation")
                    has_mesylate_to_thioether = True

                # Check for thioether to sulfoxide transformation
                if (
                    any(
                        mol.HasSubstructMatch(thioether_pattern)
                        for mol in reactant_mols
                    )
                    and not any(
                        mol.HasSubstructMatch(sulfoxide_pattern)
                        for mol in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(sulfoxide_pattern)
                ):
                    print("Detected thioether to sulfoxide transformation")
                    has_thioether_to_sulfoxide = True

                # Check for sulfoxide to sulfone transformation
                if (
                    any(
                        mol.HasSubstructMatch(sulfoxide_pattern)
                        for mol in reactant_mols
                    )
                    and not any(
                        mol.HasSubstructMatch(sulfone_pattern) for mol in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(sulfone_pattern)
                ):
                    print("Detected sulfoxide to sulfone transformation")
                    has_sulfoxide_to_sulfone = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the complete strategy is present
    strategy_present = (
        has_alcohol_to_mesylate
        and has_mesylate_to_thioether
        and has_thioether_to_sulfoxide
        and has_sulfoxide_to_sulfone
    )

    print(f"Sequential sulfur oxidation strategy detected: {strategy_present}")
    return strategy_present
