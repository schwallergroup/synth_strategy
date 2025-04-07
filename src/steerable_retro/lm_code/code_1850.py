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
    Detects a strategy involving multiple different phenol functionalization steps.
    """
    phenol_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Define SMARTS for phenol and various protected forms
                phenol_pattern = Chem.MolFromSmarts("[cH0]([OH])")
                benzyl_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][c]1[cH][cH][cH][cH][cH]1")
                mesylate_pattern = Chem.MolFromSmarts("[c][O][S](=[O])(=[O])[CH3]")
                alkyl_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][C]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and any(reactant_mols):
                    # Check for phenol in reactants
                    has_phenol = any(
                        r and r.HasSubstructMatch(phenol_pattern) for r in reactant_mols if r
                    )

                    # Check for protected forms in product
                    has_benzyl = product_mol.HasSubstructMatch(benzyl_ether_pattern)
                    has_mesylate = product_mol.HasSubstructMatch(mesylate_pattern)
                    has_alkyl_ether = product_mol.HasSubstructMatch(alkyl_ether_pattern)

                    if has_phenol and (has_benzyl or has_mesylate or has_alkyl_ether):
                        reaction_type = (
                            "benzyl_protection"
                            if has_benzyl
                            else "mesylation"
                            if has_mesylate
                            else "alkylation"
                        )
                        phenol_reactions.append((depth, reaction_type))
                        print(f"Found phenol {reaction_type} at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have at least 2 different types of phenol functionalization
    reaction_types = set(r_type for _, r_type in phenol_reactions)
    if len(reaction_types) >= 2:
        print(f"Found multiple phenol functionalization types: {reaction_types}")
        return True

    return False
