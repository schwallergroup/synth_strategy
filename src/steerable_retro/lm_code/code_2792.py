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
    This function detects a specific functional group transformation sequence:
    nitro group → amine → amide in the synthetic route.
    """
    # Track reactions with specific transformations
    nitro_to_amine_reactions = []
    amine_to_amide_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol:
                    # Patterns
                    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                    # Check for nitro to amine transformation
                    reactant_has_nitro = any(
                        mol and mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols
                    )
                    product_has_amine = product_mol.HasSubstructMatch(amine_pattern)

                    if reactant_has_nitro and product_has_amine:
                        nitro_to_amine_reactions.append((depth, node))

                    # Check for amine to amide transformation
                    reactant_has_amine = any(
                        mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
                    )

                    # Count amide bonds in reactants and product
                    reactant_amide_count = sum(
                        [
                            len(mol.GetSubstructMatches(amide_pattern)) if mol else 0
                            for mol in reactant_mols
                        ]
                    )
                    product_amide_count = (
                        len(product_mol.GetSubstructMatches(amide_pattern)) if product_mol else 0
                    )

                    if reactant_has_amine and product_amide_count > reactant_amide_count:
                        amine_to_amide_reactions.append((depth, node))

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both transformations and they're in the correct order
    # (nitro→amine should happen before amine→amide in the synthetic direction,
    # which means higher depth in retrosynthetic direction)
    if nitro_to_amine_reactions and amine_to_amide_reactions:
        # Get the depths
        nitro_to_amine_depths = [d for d, _ in nitro_to_amine_reactions]
        amine_to_amide_depths = [d for d, _ in amine_to_amide_reactions]

        # Check if any nitro→amine happens at higher depth than any amine→amide
        if any(
            n_depth > a_depth
            for n_depth in nitro_to_amine_depths
            for a_depth in amine_to_amide_depths
        ):
            print("Found nitro→amine→amide transformation sequence")
            return True

    return False
