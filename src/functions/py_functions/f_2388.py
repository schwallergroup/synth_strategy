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
    Detects if the synthesis involves aromatic iodination as a key step for enabling coupling reactions.
    """
    iodination_found = False

    def dfs_traverse(node, depth=0):
        nonlocal iodination_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an iodination reaction
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    aryl_iodide_pattern = Chem.MolFromSmarts("c-[#53]")
                    if product_mol.HasSubstructMatch(aryl_iodide_pattern):
                        # Check if reactants don't have iodine
                        has_iodine_source = False
                        has_non_iodinated_aryl = False

                        for reactant in reactants:
                            r_mol = Chem.MolFromSmiles(reactant)
                            if r_mol:
                                if (
                                    "I" in reactant
                                    or "[I]" in reactant
                                    or "Cl[I]" in reactant
                                ):
                                    has_iodine_source = True

                                if r_mol.HasSubstructMatch(
                                    Chem.MolFromSmarts("c")
                                ) and not r_mol.HasSubstructMatch(aryl_iodide_pattern):
                                    has_non_iodinated_aryl = True

                        if has_iodine_source and has_non_iodinated_aryl:
                            print("Aromatic iodination detected at depth", depth)
                            iodination_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return iodination_found
