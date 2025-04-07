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
    Detects if the synthesis route involves a convergent approach with multiple fragment couplings.
    """
    fragment_couplings = 0

    def dfs_traverse(node, depth=0):
        nonlocal fragment_couplings

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If there are multiple reactants, it might be a fragment coupling
            if len(reactants_smiles) >= 2:
                # Check if both reactants are complex (not just simple reagents)
                complex_reactants = 0
                for r_smiles in reactants_smiles:
                    try:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if (
                            r_mol and r_mol.GetNumAtoms() > 6
                        ):  # Arbitrary threshold for "complex"
                            complex_reactants += 1
                    except:
                        pass

                if complex_reactants >= 2:
                    fragment_couplings += 1
                    print(f"Fragment coupling detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return (
        fragment_couplings >= 2
    )  # Return True if at least 2 fragment couplings are detected
