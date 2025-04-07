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
    This function detects the formation of a diaryl ether (c-O-c) bond
    in the synthetic route.
    """
    has_diaryl_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_diaryl_ether_formation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for diaryl ether formation
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol:
                    # Diaryl ether pattern
                    diaryl_ether_pattern = Chem.MolFromSmarts("c[O]c")

                    # Check if product has diaryl ether
                    if product_mol.HasSubstructMatch(diaryl_ether_pattern):
                        # Check if any reactant has the same pattern
                        reactant_has_pattern = any(
                            mol and mol.HasSubstructMatch(diaryl_ether_pattern)
                            for mol in reactant_mols
                        )

                        # If product has pattern but reactants don't, it was formed in this reaction
                        if not reactant_has_pattern:
                            print(f"Found diaryl ether formation at depth {depth}")
                            has_diaryl_ether_formation = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_diaryl_ether_formation
