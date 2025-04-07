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
    This function detects if the synthesis employs a late-stage N-alkylation
    of a pyrazole ring as the final or penultimate step.
    """
    n_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction" and depth <= 1:  # Check only final or penultimate steps
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an N-alkylation of pyrazole
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol and prod_mol:
                            # Pyrazole pattern
                            pyrazole_patt = Chem.MolFromSmarts("[n]1[n]ccc1")

                            # Check if reactant has pyrazole with NH
                            nh_pyrazole_patt = Chem.MolFromSmarts("[nH]1[n]ccc1")
                            # Check if product has N-alkylated pyrazole
                            n_alkyl_pyrazole_patt = Chem.MolFromSmarts("[n]1([C])[n]ccc1")

                            if react_mol.HasSubstructMatch(
                                nh_pyrazole_patt
                            ) and prod_mol.HasSubstructMatch(n_alkyl_pyrazole_patt):
                                print(f"Late-stage pyrazole N-alkylation detected at depth {depth}")
                                n_alkylation_detected = True
                except Exception as e:
                    print(f"Error in SMILES processing: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return n_alkylation_detected
