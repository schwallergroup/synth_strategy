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
    This function detects if the synthetic route employs a thioether (C-S-C) formation
    as a key fragment coupling strategy.
    """
    thioether_formation_found = False

    def dfs_traverse(node):
        nonlocal thioether_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a coupling reaction (multiple reactants)
            if len(reactants) >= 2:
                # Check for thiol and benzyl chloride reactants
                thiol_pattern = Chem.MolFromSmarts("[SH]-[#6]")
                benzyl_cl_pattern = Chem.MolFromSmarts("[#6]-[CH2]-[Cl]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check if reactants contain thiol and benzyl chloride
                has_thiol = any(
                    mol is not None and mol.HasSubstructMatch(thiol_pattern)
                    for mol in reactant_mols
                )
                has_benzyl_cl = any(
                    mol is not None and mol.HasSubstructMatch(benzyl_cl_pattern)
                    for mol in reactant_mols
                )

                # Check if product contains thioether
                thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")
                has_thioether = product_mol is not None and product_mol.HasSubstructMatch(
                    thioether_pattern
                )

                if has_thiol and has_benzyl_cl and has_thioether:
                    print("Found thioether formation via coupling of thiol and benzyl chloride")
                    thioether_formation_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return thioether_formation_found
