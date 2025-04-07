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
    Detects a convergent synthesis strategy with Sonogashira coupling as the key fragment connection.
    Looks for C≡C bond formation between two aromatic fragments.
    """
    sonogashira_found = False
    convergent_structure = False

    def dfs_traverse(node, depth=0):
        nonlocal sonogashira_found, convergent_structure

        if node["type"] == "reaction":
            # Check if this is a Sonogashira coupling reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for C≡C bond formation
                if len(reactants) >= 2:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Look for C≡C bond in product
                        triple_bond_pattern = Chem.MolFromSmarts("[#6]#[#6]")
                        if product_mol.HasSubstructMatch(triple_bond_pattern):
                            # Check if reactants have aromatic rings and one has iodine
                            has_aryl_iodide = False
                            has_terminal_alkyne = False

                            for reactant in reactants:
                                r_mol = Chem.MolFromSmiles(reactant)
                                if r_mol:
                                    aryl_iodide_pattern = Chem.MolFromSmarts("c-[#53]")
                                    if r_mol.HasSubstructMatch(aryl_iodide_pattern):
                                        has_aryl_iodide = True

                                    terminal_alkyne_pattern = Chem.MolFromSmarts("[#6]#[#6]")
                                    if r_mol.HasSubstructMatch(terminal_alkyne_pattern):
                                        has_terminal_alkyne = True

                            if has_aryl_iodide and has_terminal_alkyne:
                                print("Sonogashira coupling detected at depth", depth)
                                sonogashira_found = True

            # Check for convergent structure - multiple children at depth > 1
            if depth > 0 and len(node.get("children", [])) > 1:
                print("Convergent structure detected at depth", depth)
                convergent_structure = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return sonogashira_found and convergent_structure
