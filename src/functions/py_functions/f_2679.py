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
    This function detects late-stage N-alkylation where two complex fragments are joined.
    """
    late_stage_n_alkylation = False

    def dfs_traverse(node):
        nonlocal late_stage_n_alkylation

        # Check only reactions at depth 0 or 1 (late stage)
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            depth = node["metadata"].get("depth", -1)
            if depth <= 1:  # Late stage reaction
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have at least 2 complex reactants
                if len(reactants) >= 2:
                    product_mol = Chem.MolFromSmiles(product)

                    # Check for N-alkylation pattern
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and product_mol:
                            # Check for nitrogen in amide/lactam
                            n_amide_pattern = Chem.MolFromSmarts("[NH][C]=O")
                            if reactant_mol.HasSubstructMatch(n_amide_pattern):
                                # Check for alkyl halide in other reactant
                                for other_reactant in reactants:
                                    if other_reactant != reactant:
                                        other_mol = Chem.MolFromSmiles(other_reactant)
                                        if other_mol:
                                            alkyl_halide = Chem.MolFromSmarts(
                                                "[C][Br,Cl,I]"
                                            )
                                            if other_mol.HasSubstructMatch(
                                                alkyl_halide
                                            ):
                                                late_stage_n_alkylation = True
                                                print(
                                                    "Late-stage N-alkylation detected"
                                                )
                                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return late_stage_n_alkylation
