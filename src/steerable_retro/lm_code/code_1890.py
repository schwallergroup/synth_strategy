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
    Detects if the synthetic route involves a late-stage coupling of complex fragments.
    Specifically looks for amide formation in the first or second reaction.
    """
    late_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_coupling_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}) and depth <= 1:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an amide coupling reaction
            if len(reactants) >= 2:  # Need at least two reactants for coupling
                product_mol = Chem.MolFromSmiles(product)
                amide_pattern = Chem.MolFromSmarts("[#7][#6](=[#8])")

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    # Check if reactants are complex (have more than 15 atoms)
                    complex_reactants = 0
                    for r in reactants:
                        mol = Chem.MolFromSmiles(r)
                        if mol and mol.GetNumAtoms() > 15:
                            complex_reactants += 1

                    if complex_reactants >= 1:
                        print(f"Found late-stage coupling at depth {depth}")
                        late_coupling_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_coupling_found
