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
    Detects if the synthetic route involves a convergent approach with 3+ fragments.
    """
    fragment_count = 0
    is_convergent = False

    def count_fragments(node):
        nonlocal fragment_count

        if node["type"] == "mol" and node.get("in_stock", False):
            fragment_count += 1
            return

        for child in node.get("children", []):
            count_fragments(child)

    def check_convergent_steps(node):
        nonlocal is_convergent

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If a step combines 2+ complex fragments, it's likely convergent
            complex_reactants = 0
            for reactant in reactants:
                if reactant and len(reactant) > 15:  # Simple heuristic for "complex" fragments
                    complex_reactants += 1

            if complex_reactants >= 2:
                print(f"Convergent step detected with {complex_reactants} complex reactants")
                is_convergent = True

        for child in node.get("children", []):
            check_convergent_steps(child)

    count_fragments(route)
    check_convergent_steps(route)

    print(f"Total fragments in synthesis: {fragment_count}")
    return is_convergent and fragment_count >= 3
