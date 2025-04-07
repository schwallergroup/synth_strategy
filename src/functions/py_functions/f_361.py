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
    Detects a strategy involving the use of protected nitrogen heterocycles
    (particularly Boc-protected piperazine) in fragment coupling.
    """
    found_boc_piperazine = False
    found_coupling_with_protected_fragment = False

    def dfs_traverse(node, depth=0):
        nonlocal found_boc_piperazine, found_coupling_with_protected_fragment

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc-protected piperazine
                boc_piperazine_pattern = Chem.MolFromSmarts("N1CCN(C(=O)OC(C)(C)C)CC1")

                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.HasSubstructMatch(boc_piperazine_pattern):
                        found_boc_piperazine = True

                        # Check if this is a coupling reaction (multiple reactants)
                        if len(reactants) >= 2:
                            # Check if other reactant has aromatic halide (typical for coupling)
                            for other_r in reactants:
                                if other_r != r:
                                    other_mol = Chem.MolFromSmiles(other_r)
                                    if other_mol and other_mol.HasSubstructMatch(
                                        Chem.MolFromSmarts("c[F,Cl,Br,I]")
                                    ):
                                        print(
                                            f"Found coupling with Boc-protected piperazine at depth {depth}"
                                        )
                                        found_coupling_with_protected_fragment = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if strategy criteria are met
    return found_boc_piperazine and found_coupling_with_protected_fragment
