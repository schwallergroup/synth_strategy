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
    Detects if the synthesis involves disconnection of a C-N bond between a piperidine ring and an amine
    """
    found_piperidine_cn_disconnection = False

    def dfs_traverse(node):
        nonlocal found_piperidine_cn_disconnection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for piperidine pattern in one fragment
            piperidine_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#7][#6][#6]1")
            amine_pattern = Chem.MolFromSmarts("[#7;H1]")

            # Check if product has both piperidine and amine connected
            product_mol = Chem.MolFromSmiles(product)

            if (
                product_mol
                and product_mol.HasSubstructMatch(piperidine_pattern)
                and product_mol.HasSubstructMatch(amine_pattern)
            ):
                # Check if reactants are split into piperidine and amine fragments
                piperidine_found = False
                amine_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(piperidine_pattern):
                            piperidine_found = True
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            amine_found = True

                if piperidine_found and amine_found:
                    print("Found piperidine-amine disconnection")
                    found_piperidine_cn_disconnection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_piperidine_cn_disconnection
