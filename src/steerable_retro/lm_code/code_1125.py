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
    This function detects if the synthesis includes multiple nucleophilic aromatic substitutions.
    """
    nas_reactions = 0

    def dfs_traverse(node):
        nonlocal nas_reactions

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for chlorinated aromatic in reactants
                chloro_aromatic_found = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        chloro_aromatic_pattern = Chem.MolFromSmarts("[c][Cl]")
                        if reactant_mol.HasSubstructMatch(chloro_aromatic_pattern):
                            chloro_aromatic_found = True
                            break

                # Check for new C-N bond in product where chlorine was
                if chloro_aromatic_found:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # This is a simplified check - in reality would need to confirm
                        # the exact position where Cl was replaced
                        c_n_bond_pattern = Chem.MolFromSmarts("[c][N]")
                        if product_mol.HasSubstructMatch(c_n_bond_pattern):
                            nas_reactions += 1
                            print(
                                f"Detected nucleophilic aromatic substitution, count: {nas_reactions}"
                            )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nas_reactions >= 2  # Return True if at least 2 NAS reactions are found
