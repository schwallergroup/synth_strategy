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
    This function detects a strategy involving heterocycle formation (imidazole and pyrazole)
    combined with nitro group chemistry (introduction and reduction).
    """
    # Track if we found the key elements of the strategy
    has_nitro_introduction = False
    has_nitro_reduction = False
    has_imidazole_formation = False
    has_pyrazole_formation = False

    def dfs_traverse(node):
        nonlocal has_nitro_introduction, has_nitro_reduction, has_imidazole_formation, has_pyrazole_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Convert to RDKit molecules
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product and all(r for r in reactants):
                        # Check for nitro introduction
                        nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
                        reactant_nitro_count = sum(
                            len(r.GetSubstructMatches(nitro_pattern))
                            for r in reactants
                            if r
                        )
                        product_nitro_count = len(
                            product.GetSubstructMatches(nitro_pattern)
                        )

                        if product_nitro_count > reactant_nitro_count:
                            has_nitro_introduction = True
                            print("Detected nitro group introduction")

                        # Check for nitro reduction
                        amine_pattern = Chem.MolFromSmarts("[NH2]")
                        reactant_amine_count = sum(
                            len(r.GetSubstructMatches(amine_pattern))
                            for r in reactants
                            if r
                        )
                        product_amine_count = len(
                            product.GetSubstructMatches(amine_pattern)
                        )

                        if product_amine_count > reactant_amine_count and any(
                            len(r.GetSubstructMatches(nitro_pattern)) > 0
                            for r in reactants
                            if r
                        ):
                            has_nitro_reduction = True
                            print("Detected nitro reduction to amine")

                        # Check for imidazole formation
                        imidazole_pattern = Chem.MolFromSmarts("[nH]1cncc1")
                        reactant_imidazole_count = sum(
                            len(r.GetSubstructMatches(imidazole_pattern))
                            for r in reactants
                            if r
                        )
                        product_imidazole_count = len(
                            product.GetSubstructMatches(imidazole_pattern)
                        )

                        if product_imidazole_count > reactant_imidazole_count:
                            has_imidazole_formation = True
                            print("Detected imidazole formation")

                        # Check for pyrazole formation
                        pyrazole_pattern = Chem.MolFromSmarts("[nH]1ncc[c]1")
                        reactant_pyrazole_count = sum(
                            len(r.GetSubstructMatches(pyrazole_pattern))
                            for r in reactants
                            if r
                        )
                        product_pyrazole_count = len(
                            product.GetSubstructMatches(pyrazole_pattern)
                        )

                        if product_pyrazole_count > reactant_pyrazole_count:
                            has_pyrazole_formation = True
                            print("Detected pyrazole formation")

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # The strategy is present if we have heterocycle formation and nitro chemistry
    strategy_present = (has_imidazole_formation or has_pyrazole_formation) and (
        has_nitro_introduction or has_nitro_reduction
    )

    if strategy_present:
        print("Detected heterocycle formation with nitro chemistry strategy")

    return strategy_present
