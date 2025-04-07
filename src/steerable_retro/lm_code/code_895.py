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
    This function detects a late-stage Suzuki coupling strategy (depth 0-1)
    involving a halogenated aryl partner.
    """
    suzuki_coupling_found = False
    halogenated_partner = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_found, halogenated_partner

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0-1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid derivative in reactants
                boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")

                # Check for aryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("c-[#35,#53,#17]")

                # Check for fluorinated aryl group
                fluorinated_aryl = Chem.MolFromSmarts("c[F]")

                # Check if this looks like a Suzuki coupling
                has_boronic = False
                has_aryl_halide = False
                has_fluorinated_aryl = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(boronic_pattern):
                            has_boronic = True
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                        if mol.HasSubstructMatch(fluorinated_aryl):
                            has_fluorinated_aryl = True

                # Product should have a new C-C bond
                if has_boronic and has_aryl_halide:
                    suzuki_coupling_found = True
                    if has_fluorinated_aryl:
                        halogenated_partner = True
                    print(f"Found Suzuki coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found a late-stage Suzuki coupling with halogenated partner
    return suzuki_coupling_found and halogenated_partner
