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
    This function detects a protection-deprotection sequence using trifluoroacetyl group.
    It looks for protection of an amine with trifluoroacetyl group and subsequent deprotection.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for trifluoroacetyl protection (amine to trifluoroacetamide)
            if not protection_found:
                amine_pattern = Chem.MolFromSmarts("[NH2][c]")
                trifluoroacetamide_pattern = Chem.MolFromSmarts("[NH]C(=O)C(F)(F)(F)")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        trifluoroacetamide_pattern
                    ):
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                amine_pattern
                            ):
                                protection_found = True
                                print("Trifluoroacetyl protection detected")
                                break
                except:
                    pass

            # Check for trifluoroacetyl deprotection (trifluoroacetamide to amine)
            if not deprotection_found:
                trifluoroacetamide_pattern = Chem.MolFromSmarts("[NH]C(=O)C(F)(F)(F)")
                amine_pattern = Chem.MolFromSmarts("[NH2][c]")

                try:
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            trifluoroacetamide_pattern
                        ):
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and product_mol.HasSubstructMatch(
                                amine_pattern
                            ):
                                deprotection_found = True
                                print("Trifluoroacetyl deprotection detected")
                                break
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return protection_found and deprotection_found
