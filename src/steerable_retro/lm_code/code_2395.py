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
    This function detects synthetic routes involving polyethylene glycol (PEG) chain
    manipulation, including chain cleavage and terminal functional group conversion.
    """
    # Track if we found evidence of PEG manipulation
    has_polyether_chain = False
    has_terminal_group_conversion = False

    def dfs_traverse(node):
        nonlocal has_polyether_chain, has_terminal_group_conversion

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            products = parts[2].split(".")

            # Check for polyether chains in reactants or products
            for mol_smiles in reactants + products:
                try:
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        # Check for polyether pattern (at least 2 consecutive ether linkages)
                        peg_pattern = Chem.MolFromSmarts("[#8]-[#6]-[#6]-[#8]-[#6]-[#6]")
                        if mol.HasSubstructMatch(peg_pattern):
                            has_polyether_chain = True
                            print(f"Found polyether chain in: {mol_smiles}")
                except:
                    continue

            # Check for terminal group conversion (alcohol to halide)
            for reactant in reactants:
                try:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol:
                        # Check for terminal alcohol
                        terminal_alcohol = Chem.MolFromSmarts("[#6]-[#8;H1]")
                        if r_mol.HasSubstructMatch(terminal_alcohol):
                            # Look for corresponding product with terminal halide
                            for product in products:
                                try:
                                    p_mol = Chem.MolFromSmiles(product)
                                    if p_mol:
                                        terminal_halide = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")
                                        if p_mol.HasSubstructMatch(terminal_halide):
                                            has_terminal_group_conversion = True
                                            print(
                                                f"Found terminal group conversion: {reactant} -> {product}"
                                            )
                                except:
                                    continue
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found evidence of PEG manipulation strategy
    strategy_detected = has_polyether_chain and has_terminal_group_conversion
    print(f"Polyether chain manipulation strategy detected: {strategy_detected}")
    return strategy_detected
