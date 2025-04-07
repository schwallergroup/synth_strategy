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
    Detects a Boc protection/deprotection sequence in the synthesis.
    """
    # Track if we found Boc protection and deprotection
    boc_protected_found = False
    boc_deprotection_found = False

    def dfs_traverse(node):
        nonlocal boc_protected_found, boc_deprotection_found

        if node["type"] == "mol":
            # Check if molecule contains Boc-protected amine
            if node.get("smiles"):
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(
                    Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
                ):
                    boc_protected_found = True
                    print("Found Boc-protected intermediate")

        elif node["type"] == "reaction":
            # Check if this is a Boc deprotection reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant has Boc but product doesn't
                reactant_has_boc = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(
                        Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
                    ):
                        reactant_has_boc = True

                product_mol = Chem.MolFromSmiles(product)
                product_has_boc = product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
                )

                if reactant_has_boc and not product_has_boc:
                    boc_deprotection_found = True
                    print("Found Boc deprotection step")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return boc_protected_found and boc_deprotection_found
