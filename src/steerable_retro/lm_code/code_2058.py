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
    This function detects a strategy involving transformation of a carbonyl group to a chloride,
    which is then used for subsequent substitution reactions.
    """
    carbonyl_to_chloride = False
    chloride_substitution = False

    def dfs_traverse(node, depth=0):
        nonlocal carbonyl_to_chloride, chloride_substitution

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and all(r is not None for r in reactant_mols):
                # Check for carbonyl to chloride transformation
                carbonyl_pattern = Chem.MolFromSmarts("[c](=[O])")
                chloride_pattern = Chem.MolFromSmarts("[c]-[Cl]")

                reactants_have_carbonyl = any(
                    r.HasSubstructMatch(carbonyl_pattern) for r in reactant_mols if r
                )
                product_has_chloride = (
                    product_mol.HasSubstructMatch(chloride_pattern) if product_mol else False
                )

                if reactants_have_carbonyl and product_has_chloride:
                    print(f"Carbonyl to chloride transformation detected at depth {depth}")
                    carbonyl_to_chloride = True

                # Check for chloride substitution
                amine_pattern = Chem.MolFromSmarts("[c]-[NH]")
                reactants_have_chloride = any(
                    r.HasSubstructMatch(chloride_pattern) for r in reactant_mols if r
                )
                product_has_amine = (
                    product_mol.HasSubstructMatch(amine_pattern) if product_mol else False
                )

                if reactants_have_chloride and product_has_amine:
                    print(f"Chloride substitution detected at depth {depth}")
                    chloride_substitution = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if the strategy is detected
    strategy_detected = carbonyl_to_chloride and chloride_substitution
    print(
        f"Strategy detected: {strategy_detected} (Carbonyl to chloride: {carbonyl_to_chloride}, Chloride substitution: {chloride_substitution})"
    )
    return strategy_detected
