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
    This function detects a synthetic strategy involving the formation of a chroman ring system
    through intramolecular cyclization.
    """
    chroman_formed = False

    def dfs_traverse(node):
        nonlocal chroman_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Check if chroman pattern exists in product but not in reactants
                    chroman_pattern = Chem.MolFromSmarts(
                        "[c]1[c][c][c][c][c]1[O][C][C][C]"
                    )
                    if product_mol and chroman_pattern:
                        if product_mol.HasSubstructMatch(chroman_pattern):
                            # Check if pattern doesn't exist in any reactant
                            if not any(
                                r and r.HasSubstructMatch(chroman_pattern)
                                for r in reactant_mols
                                if r
                            ):
                                chroman_formed = True
                                print("Chroman ring formation detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return chroman_formed
