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
    Detects if the synthetic route involves a Suzuki coupling (aryl-aryl C-C bond formation
    between an aryl bromide/iodide and a boronate ester).
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant contains boronate and another contains Br/I on aromatic
            has_boronate = False
            has_aryl_halide = False

            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for boronate
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[B]")):
                            has_boronate = True
                        # Check for aryl bromide/iodide
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[Br,I]")):
                            has_aryl_halide = True

            # Check if product has biaryl bond that wasn't in reactants
            if has_boronate and has_aryl_halide:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[c]-[c]")
                ):
                    print("Suzuki coupling detected")
                    suzuki_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected
