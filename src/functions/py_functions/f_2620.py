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
    This function detects an amide formation strategy using an acid chloride intermediate.
    """
    # Track if we found the strategy components
    found_acid_chloride = False
    found_amide_formation = False

    # SMARTS patterns
    acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")
    amine_pattern = Chem.MolFromSmarts("[NH]")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")

    def dfs_traverse(node, depth=0):
        nonlocal found_acid_chloride, found_amide_formation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    product_mol = Chem.MolFromSmiles(product)

                    # Check for acid chloride formation
                    if any(
                        mol is not None and mol.HasSubstructMatch(acid_chloride_pattern)
                        for mol in reactant_mols
                    ):
                        found_acid_chloride = True
                        print(f"Found acid chloride at depth {depth}")

                    # Check for amide formation from acid chloride
                    has_acid_chloride = any(
                        mol is not None and mol.HasSubstructMatch(acid_chloride_pattern)
                        for mol in reactant_mols
                    )
                    has_amine = any(
                        mol is not None and mol.HasSubstructMatch(amine_pattern)
                        for mol in reactant_mols
                    )
                    has_amide_product = (
                        product_mol is not None
                        and product_mol.HasSubstructMatch(amide_pattern)
                    )

                    if has_acid_chloride and has_amine and has_amide_product:
                        found_amide_formation = True
                        print(f"Found amide formation at depth {depth}")
                except:
                    print("Error processing SMILES in reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we found both acid chloride and amide formation
    strategy_present = found_acid_chloride and found_amide_formation
    print(f"Amide formation from acid chloride: {strategy_present}")
    return strategy_present
