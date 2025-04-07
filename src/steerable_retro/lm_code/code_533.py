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
    Detects if the route uses a protection-deprotection sequence,
    specifically looking for N-acetylation followed by deacetylation.
    """
    acetylation_depths = []
    deacetylation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal acetylation_depths, deacetylation_depths

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acetylation: NH2 + something → NHC(=O)CH3
            acetylation_pattern = Chem.MolFromSmarts("NC(=O)C")
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(acetylation_pattern):
                # Check if pattern is not in reactants
                pattern_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(acetylation_pattern):
                        pattern_in_reactants = True
                        break

                if not pattern_in_reactants:
                    acetylation_depths.append(depth)
                    print(f"Detected acetylation at depth {depth}")

            # Check for deacetylation: NHC(=O)CH3 → NH2
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                # Check if acetyl group was in reactants
                acetyl_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("NC(=O)C")
                    ):
                        acetyl_in_reactants = True
                        break

                if acetyl_in_reactants:
                    deacetylation_depths.append(depth)
                    print(f"Detected deacetylation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there's at least one acetylation followed by deacetylation
    for acetyl_depth in acetylation_depths:
        for deacetyl_depth in deacetylation_depths:
            if deacetyl_depth < acetyl_depth:  # Remember: lower depth = later in synthesis
                print(
                    f"Protection-deprotection sequence detected: acetylation at depth {acetyl_depth}, deacetylation at depth {deacetyl_depth}"
                )
                return True

    return False
