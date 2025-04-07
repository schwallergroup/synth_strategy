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
    This function detects a synthetic strategy involving THP protection of an alcohol
    followed by its deprotection later in the synthesis.
    """
    # Track if we've seen THP protection and deprotection
    protection_found = False
    deprotection_found = False

    # SMARTS for THP protected alcohol
    thp_pattern = Chem.MolFromSmarts("[#8]-[CH]1[CH2][CH2][CH2][CH2][O]1")

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                # Try to create RDKit molecules
                try:
                    reactants = [
                        Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")
                    ]
                    product = Chem.MolFromSmiles(products_smiles)

                    # Check for THP protection reaction
                    if product and any(r for r in reactants if r):
                        # Protection: If THP appears in product but not in reactants
                        reactants_have_thp = any(
                            r.HasSubstructMatch(thp_pattern) for r in reactants if r
                        )
                        product_has_thp = product.HasSubstructMatch(thp_pattern)

                        if product_has_thp and not reactants_have_thp:
                            # Check if one reactant has dihydropyran pattern
                            dihydropyran_pattern = Chem.MolFromSmarts(
                                "[CH]1=[CH][CH2][CH2][CH2][O]1"
                            )
                            if any(
                                r.HasSubstructMatch(dihydropyran_pattern)
                                for r in reactants
                                if r
                            ):
                                protection_found = True
                                print("THP protection detected")

                        # Deprotection: If THP appears in reactants but not in product
                        if reactants_have_thp and not product_has_thp:
                            deprotection_found = True
                            print("THP deprotection detected")
                except:
                    pass  # Skip if SMILES parsing fails

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return protection_found and deprotection_found
