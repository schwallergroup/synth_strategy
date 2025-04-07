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
    Detects if the synthetic route contains SNAr reactions on fluorinated aromatics.
    """
    snar_count = 0

    def dfs_traverse(node):
        nonlocal snar_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an SNAr reaction on fluorinated aromatics
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Check for fluorinated aromatic in reactant
                    fluoro_aromatic_pattern = Chem.MolFromSmarts("c[F]")
                    if reactant_mol.HasSubstructMatch(fluoro_aromatic_pattern):
                        # Check for nucleophile attack patterns (C-N or C-O bond formation)
                        # where a C-F bond was broken
                        if len(
                            reactant_mol.GetSubstructMatches(fluoro_aromatic_pattern)
                        ) > len(
                            product_mol.GetSubstructMatches(fluoro_aromatic_pattern)
                        ):

                            # Check for new C-N or C-O bond in product
                            c_n_pattern = Chem.MolFromSmarts("c[NH]c")  # Diarylamine
                            c_o_pattern = Chem.MolFromSmarts("c[O][CH3]")  # Methoxy

                            if product_mol.HasSubstructMatch(
                                c_n_pattern
                            ) or product_mol.HasSubstructMatch(c_o_pattern):
                                snar_count += 1
                                print(
                                    f"Detected SNAr reaction at depth {node.get('depth', 'unknown')}"
                                )
                                break

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if SNAr reactions on fluorinated aromatics were detected
    result = snar_count > 0
    print(
        f"SNAr on fluorinated aromatics strategy detected: {result} (count: {snar_count})"
    )
    return result
