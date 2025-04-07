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
    Detects if the synthesis uses a Boc protection/deprotection sequence.
    Looks for Boc group (tert-butyloxycarbonyl) in earlier steps and its removal in later steps.
    """
    boc_protected_steps = []
    boc_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protected_steps, boc_deprotection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Boc group pattern
                boc_pattern = Chem.MolFromSmarts("[#6]OC(=O)[#7]")

                # Check for Boc in reactants and product
                reactants_have_boc = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(boc_pattern):
                        reactants_have_boc = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                product_has_boc = product_mol and product_mol.HasSubstructMatch(
                    boc_pattern
                )

                # If reactants have Boc but product doesn't, it's a deprotection
                if reactants_have_boc and not product_has_boc:
                    boc_deprotection_found = True
                    print(f"Detected Boc deprotection at depth {depth}")

                # If product has Boc, record the depth
                if product_has_boc:
                    boc_protected_steps.append(depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we found both protection and deprotection
    has_sequence = boc_deprotection_found and len(boc_protected_steps) > 0
    if has_sequence:
        print(
            f"Detected Boc protection/deprotection sequence. Protected at depths: {boc_protected_steps}"
        )

    return has_sequence
