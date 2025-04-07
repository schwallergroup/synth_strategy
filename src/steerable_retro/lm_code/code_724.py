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
    This function detects if the synthetic route involves late-stage reductive amination
    (formation of secondary amine from primary amine and aldehyde).
    """
    reductive_amination_found = False
    reductive_amination_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found, reductive_amination_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols if r):
                # Check for reductive amination pattern
                aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")
                primary_amine_pattern = Chem.MolFromSmarts("[NH2]")
                secondary_amine_pattern = Chem.MolFromSmarts("[NH1][CH2]")

                reactants_have_aldehyde = any(
                    r.HasSubstructMatch(aldehyde_pattern) for r in reactant_mols if r
                )
                reactants_have_primary_amine = any(
                    r.HasSubstructMatch(primary_amine_pattern) for r in reactant_mols if r
                )
                product_has_secondary_amine = product_mol.HasSubstructMatch(secondary_amine_pattern)

                if (
                    reactants_have_aldehyde
                    and reactants_have_primary_amine
                    and product_has_secondary_amine
                ):
                    print(f"Reductive amination detected at depth {depth}")
                    reductive_amination_found = True
                    reductive_amination_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if reductive amination occurred in the late stage (depth <= 1)
    return reductive_amination_found and reductive_amination_depth <= 1
