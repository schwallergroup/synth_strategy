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
    This function detects convergent synthesis with Suzuki coupling.
    It looks for a synthesis where two complex fragments are joined via Suzuki coupling.
    """
    suzuki_found = False
    is_convergent = False

    def dfs_traverse(node):
        nonlocal suzuki_found, is_convergent

        if node["type"] == "reaction":
            # Check if this is a Suzuki coupling
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid/ester in reactants
                boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
                aryl_halide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#35,#53,#17]")

                has_boronic = False
                has_aryl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(boronic_pattern):
                            has_boronic = True
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True

                # If both patterns are found and product has a new biaryl bond, it's likely a Suzuki coupling
                if has_boronic and has_aryl_halide:
                    # Check if product has a biaryl bond that wasn't in reactants
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")
                        if prod_mol.HasSubstructMatch(biaryl_pattern):
                            suzuki_found = True

                            # Check if both fragments are complex (more than 10 atoms)
                            complex_fragments = 0
                            for reactant in reactants:
                                mol = Chem.MolFromSmiles(reactant)
                                if mol and mol.GetNumAtoms() > 10:
                                    complex_fragments += 1

                            if complex_fragments >= 2:
                                is_convergent = True
                                print("Found convergent Suzuki coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return suzuki_found and is_convergent
