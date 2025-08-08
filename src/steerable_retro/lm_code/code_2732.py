#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects if N-alkylation is used as a fragment coupling strategy.
    """
    n_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for N-alkylation pattern: typically involves a halide (Cl, Br) and an amine
            if (
                ("Cl" in reactants_part or "Br" in reactants_part)
                and "N" in reactants_part
                and "N" in product_part
            ):

                # More specific check: look for C-N bond formation
                reactants = reactants_part.split(".")
                if len(reactants) >= 2:  # Need at least two fragments for coupling
                    # One fragment should have halide, one should have nitrogen
                    halide_fragment = False
                    amine_fragment = False

                    for reactant in reactants:
                        if "Cl" in reactant or "Br" in reactant:
                            halide_fragment = True
                        if "N" in reactant:
                            amine_fragment = True

                    if halide_fragment and amine_fragment:
                        n_alkylation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"N-alkylation fragment coupling: {n_alkylation_detected}")
    return n_alkylation_detected
