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
    This function detects if phosphonate groups are preserved throughout the synthesis.
    """
    phosphonate_pattern = Chem.MolFromSmarts("[P](=[O])([O][#6])[O][#6]")

    # Track if phosphonates are present and preserved
    phosphonate_present = False
    phosphonate_modified = False

    def dfs_traverse(node):
        nonlocal phosphonate_present, phosphonate_modified

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Count phosphonate groups in reactants and product
                reactant_phosphonates = sum(
                    len(mol.GetSubstructMatches(phosphonate_pattern))
                    for mol in reactant_mols
                    if mol
                )

                product_phosphonates = (
                    len(product_mol.GetSubstructMatches(phosphonate_pattern))
                    if product_mol
                    else 0
                )

                if reactant_phosphonates > 0 or product_phosphonates > 0:
                    phosphonate_present = True

                # If phosphonate count decreases, they're being modified
                if product_phosphonates < reactant_phosphonates:
                    phosphonate_modified = True
                    print("Phosphonate groups are being modified")
            except:
                print("Error processing molecules for phosphonate detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if phosphonates are present and preserved
    return phosphonate_present and not phosphonate_modified
