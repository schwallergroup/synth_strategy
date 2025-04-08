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
    Detects if the synthetic route includes a nucleophilic aromatic substitution
    (displacement of halide by amine).
    """
    has_snar = False

    def dfs_traverse(node):
        nonlocal has_snar

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Look for a reaction where one reactant has an aryl halide and another has an amine
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#9,#17,#35,#53]")
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

                # Check if product has aryl amine
                aryl_amine_pattern = Chem.MolFromSmarts("[c]-[#7]")

                has_aryl_halide = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(aryl_amine_pattern):
                    if has_aryl_halide and has_amine:
                        has_snar = True
                        print(f"Found nucleophilic aromatic substitution: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_snar
