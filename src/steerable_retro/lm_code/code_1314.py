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
    This function detects N-arylation of a pyrazole heterocycle.
    """
    pyrazole_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#7][#6]1")
    n_aryl_pyrazole_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#7](-[c])[#6]1")

    strategy_detected = False

    def dfs_traverse(node):
        nonlocal strategy_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains N-arylated pyrazole
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check if any reactant has pyrazole without N-arylation
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)

                    if (
                        product_mol
                        and reactant_mol
                        and product_mol.HasSubstructMatch(n_aryl_pyrazole_pattern)
                        and reactant_mol.HasSubstructMatch(pyrazole_pattern)
                        and not reactant_mol.HasSubstructMatch(n_aryl_pyrazole_pattern)
                    ):
                        print("Detected N-arylation of pyrazole")
                        strategy_detected = True
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return strategy_detected
