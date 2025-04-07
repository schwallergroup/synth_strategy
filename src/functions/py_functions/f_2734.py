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
    This function detects a synthetic strategy involving late-stage amide coupling
    using an acyl chloride intermediate.
    """
    # Track if we found the amide coupling
    found_amide_coupling = False

    def dfs_traverse(node):
        nonlocal found_amide_coupling

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is the root reaction (depth 0)
            if node.get("depth", 0) == 0:
                # Check for amide formation from acyl chloride
                acyl_chloride_pattern = Chem.MolFromSmarts("[C$(C=O)][Cl]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                amide_pattern = Chem.MolFromSmarts("[NH][C$(C=O)]")

                # Check reactants for acyl chloride and amine
                has_acyl_chloride = False
                has_amine = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(acyl_chloride_pattern):
                        has_acyl_chloride = True
                    if mol and mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

                # Check product for amide
                product_mol = Chem.MolFromSmiles(product)
                has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

                if has_acyl_chloride and has_amine and has_amide:
                    found_amide_coupling = True
                    print("Found late-stage amide coupling using acyl chloride")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_amide_coupling
