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
    This function detects N-alkylation as a fragment coupling strategy.
    """
    n_alkylation_found = False

    def dfs_traverse(node):
        nonlocal n_alkylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for N-alkylation
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
                mesylate_pattern = Chem.MolFromSmarts("[#6]-[#8][#16](=[#8])(=[#8])[#6]")

                product_mol = Chem.MolFromSmiles(product)

                amine_found = False
                mesylate_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            amine_found = True
                        if reactant_mol.HasSubstructMatch(mesylate_pattern):
                            mesylate_found = True

                if amine_found and mesylate_found:
                    n_alkylation_found = True
                    print("N-alkylation fragment coupling detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return n_alkylation_found
