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
    This function detects if borylation of aryl halide is used for cross-coupling.
    """
    borylation_found = False

    def dfs_traverse(node):
        nonlocal borylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide to aryl boronic ester transformation
                aryl_halide_pattern = Chem.MolFromSmarts("c-[#35,#53,#17]")
                diboron_pattern = Chem.MolFromSmarts("[#5]-[#8,#5]")
                aryl_boron_pattern = Chem.MolFromSmarts("c-[#5](-[#8])-[#8]")

                has_aryl_halide = False
                has_diboron = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                        if mol.HasSubstructMatch(diboron_pattern):
                            has_diboron = True

                # Check if product has aryl-boron bond
                prod_mol = Chem.MolFromSmiles(product)
                has_aryl_boron = prod_mol and prod_mol.HasSubstructMatch(aryl_boron_pattern)

                if has_aryl_halide and has_diboron and has_aryl_boron:
                    borylation_found = True
                    print("Found borylation of aryl halide")

                    # Check if this borylation is followed by a coupling reaction
                    for child in node.get("children", []):
                        if child["type"] == "reaction":
                            if "rsmi" in child.get("metadata", {}):
                                child_rsmi = child["metadata"]["rsmi"]
                                child_reactants = child_rsmi.split(">")[0].split(".")
                                child_product = child_rsmi.split(">")[-1]

                                # Check for Suzuki coupling patterns
                                has_aryl_halide_in_child = False
                                for reactant in child_reactants:
                                    mol = Chem.MolFromSmiles(reactant)
                                    if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                                        has_aryl_halide_in_child = True

                                if has_aryl_halide_in_child and has_aryl_boron:
                                    print("Borylation is followed by cross-coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return borylation_found
