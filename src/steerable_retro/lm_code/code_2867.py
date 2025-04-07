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
    Detects if the synthetic route contains a late-stage N-methylation via reductive amination.
    Specifically looks for addition of a methyl group to a nitrogen atom in the final steps.
    """
    found_n_methylation = False

    def is_n_methylation(reactants, product):
        # Check if one reactant is formaldehyde
        formaldehyde_pattern = Chem.MolFromSmarts("[CH2]=O")
        has_formaldehyde = False
        for r in reactants:
            r_mol = Chem.MolFromSmiles(r)
            if r_mol and r_mol.HasSubstructMatch(formaldehyde_pattern):
                has_formaldehyde = True
                break

        if not has_formaldehyde:
            return False

        # Check if another reactant has NH or NH2 group
        amine_pattern = Chem.MolFromSmarts("[N;H1,H2]")
        has_amine = False
        for r in reactants:
            r_mol = Chem.MolFromSmiles(r)
            if r_mol and r_mol.HasSubstructMatch(amine_pattern):
                has_amine = True
                break

        if not has_amine:
            return False

        # Check if product has N-CH3 group
        n_methyl_pattern = Chem.MolFromSmarts("[N]-[CH3]")
        product_mol = Chem.MolFromSmiles(product)
        if product_mol and product_mol.HasSubstructMatch(n_methyl_pattern):
            return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_methylation

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check only late-stage reactions (depth 0 or 1)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            if is_n_methylation(reactants, product):
                found_n_methylation = True
                print(f"Found N-methylation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage N-methylation: {found_n_methylation}")
    return found_n_methylation
