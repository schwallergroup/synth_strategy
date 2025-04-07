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
    This function detects a synthetic strategy involving multiple C-N bond
    formations, including at least one N-arylation step.
    """
    c_n_bond_formations = 0
    n_arylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal c_n_bond_formations, n_arylation_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for C-N bond formation
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Look for new C-N bonds
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        # Compare product and reactant for new C-N bonds
                        prod_bonds = set(
                            [
                                (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                                for b in product_mol.GetBonds()
                                if (
                                    product_mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetSymbol()
                                    == "C"
                                    and product_mol.GetAtomWithIdx(b.GetEndAtomIdx()).GetSymbol()
                                    == "N"
                                )
                                or (
                                    product_mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetSymbol()
                                    == "N"
                                    and product_mol.GetAtomWithIdx(b.GetEndAtomIdx()).GetSymbol()
                                    == "C"
                                )
                            ]
                        )

                        react_bonds = set(
                            [
                                (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                                for b in reactant_mol.GetBonds()
                                if (
                                    reactant_mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetSymbol()
                                    == "C"
                                    and reactant_mol.GetAtomWithIdx(b.GetEndAtomIdx()).GetSymbol()
                                    == "N"
                                )
                                or (
                                    reactant_mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetSymbol()
                                    == "N"
                                    and reactant_mol.GetAtomWithIdx(b.GetEndAtomIdx()).GetSymbol()
                                    == "C"
                                )
                            ]
                        )

                        if len(prod_bonds) > len(react_bonds):
                            c_n_bond_formations += 1
                            print(f"C-N bond formation detected at depth {depth}")

                            # Check if it's an N-arylation (aromatic C-N bond)
                            if "c" in product and "I" in "".join(reactants):
                                n_arylation_pattern = Chem.MolFromSmarts("c-[#7]")
                                if product_mol.HasSubstructMatch(n_arylation_pattern):
                                    n_arylation_detected = True
                                    print(f"N-arylation detected at depth {depth}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have multiple C-N bond formations including at least one N-arylation
    return c_n_bond_formations >= 2 and n_arylation_detected
