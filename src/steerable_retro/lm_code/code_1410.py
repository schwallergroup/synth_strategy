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
    This function detects if the synthesis includes a transition metal-catalyzed biaryl formation.
    """
    # Track if we found biaryl coupling
    found_biaryl_coupling = False

    # SMARTS patterns
    aryl_halide_pattern = Chem.MolFromSmarts("c-[Br,I,Cl]")
    biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")  # Simplified pattern for biaryl

    def dfs_traverse(node):
        nonlocal found_biaryl_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for biaryl coupling
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                has_aryl_halide = any(
                    mol and mol.HasSubstructMatch(aryl_halide_pattern) for mol in reactant_mols
                )
                has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)

                if has_aryl_halide and has_biaryl:
                    found_biaryl_coupling = True
                    print("Found biaryl formation via coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_biaryl_coupling
