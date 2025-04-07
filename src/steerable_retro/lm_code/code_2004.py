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
    This function detects if a Friedel-Crafts acylation is present in the synthesis.
    """
    friedel_crafts_found = False

    def dfs_traverse(node):
        nonlocal friedel_crafts_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride in reactants
                acid_chloride_pattern = Chem.MolFromSmarts("[C](=O)Cl")
                # Check for aromatic in reactants
                aromatic_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1")
                # Check for diaryl ketone in product
                diaryl_ketone_pattern = Chem.MolFromSmarts("[c]-[C](=O)-[c]")

                acid_chloride_found = False
                aromatic_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(acid_chloride_pattern):
                        acid_chloride_found = True

                    if reactant_mol.HasSubstructMatch(aromatic_pattern):
                        aromatic_found = True

                product_mol = Chem.MolFromSmiles(product)
                if (
                    acid_chloride_found
                    and aromatic_found
                    and product_mol
                    and product_mol.HasSubstructMatch(diaryl_ketone_pattern)
                ):
                    friedel_crafts_found = True
                    print("Found Friedel-Crafts acylation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return friedel_crafts_found
