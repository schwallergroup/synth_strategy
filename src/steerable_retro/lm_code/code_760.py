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
    Detects the strategy of introducing a complex amine fragment (like trifluoromethyl phenyl piperazine)
    in a late-stage coupling reaction.
    """
    # Initialize tracking variables
    complex_amine_coupling = False

    # SMARTS patterns for complex amine fragments
    trifluoromethyl_pattern = Chem.MolFromSmarts("[#6][C]([F])([F])[F]")
    piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")

    def dfs_traverse(node, depth):
        nonlocal complex_amine_coupling

        if (
            node["type"] == "reaction"
            and depth <= 1
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            # Focus on late-stage reactions (low depth)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains both trifluoromethyl and piperazine
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    has_cf3 = reactant_mol.HasSubstructMatch(trifluoromethyl_pattern)
                    has_piperazine = reactant_mol.HasSubstructMatch(piperazine_pattern)

                    if has_cf3 and has_piperazine:
                        product_mol = Chem.MolFromSmiles(product)
                        if (
                            product_mol
                            and product_mol.HasSubstructMatch(trifluoromethyl_pattern)
                            and product_mol.HasSubstructMatch(piperazine_pattern)
                        ):
                            print(f"Complex amine fragment coupling detected at depth {depth}")
                            complex_amine_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route, 0)

    print(f"Late-stage fragment coupling strategy detected: {complex_amine_coupling}")
    return complex_amine_coupling
