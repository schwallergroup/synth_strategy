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

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects the presence of a fluorinated aromatic group in the final product.
    Specifically looking for difluoro or trifluoro groups.
    """
    has_fluorinated_aromatic = False

    def dfs_traverse(node, depth=0):
        nonlocal has_fluorinated_aromatic

        if node["type"] == "mol" and "smiles" in node:
            # Check if this is the final product (at root, depth 0)
            if depth == 0:
                print(f"Examining final product: {node['smiles']}")

                # Check for trifluoro group
                if checker.check_fg("Trifluoro group", node["smiles"]):
                    print(f"Found trifluoro group in final product: {node['smiles']}")
                    has_fluorinated_aromatic = True
                    return

                # Check for aromatic rings with multiple fluorine atoms
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Create a dictionary to count fluorines per aromatic ring
                    aromatic_atoms = set()
                    for atom in mol.GetAtoms():
                        if atom.GetIsAromatic():
                            aromatic_atoms.add(atom.GetIdx())

                    # Count fluorines attached to aromatic atoms
                    fluorine_count = 0
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() == "F":
                            # Check if this fluorine is attached to an aromatic atom
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetIdx() in aromatic_atoms:
                                    fluorine_count += 1
                                    break

                    if fluorine_count >= 2:
                        print(
                            f"Found multiple ({fluorine_count}) fluorines on aromatic ring in final product: {node['smiles']}"
                        )
                        has_fluorinated_aromatic = True

            # Also check all molecules in the route for debugging
            else:
                if checker.check_fg("Trifluoro group", node["smiles"]):
                    print(
                        f"Found trifluoro group in intermediate at depth {depth}: {node['smiles']}"
                    )

                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Count fluorines attached to aromatic atoms
                    aromatic_atoms = set()
                    for atom in mol.GetAtoms():
                        if atom.GetIsAromatic():
                            aromatic_atoms.add(atom.GetIdx())

                    fluorine_count = 0
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() == "F":
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetIdx() in aromatic_atoms:
                                    fluorine_count += 1
                                    break

                    if fluorine_count >= 2:
                        print(
                            f"Found multiple ({fluorine_count}) fluorines on aromatic ring at depth {depth}: {node['smiles']}"
                        )

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Fluorinated aromatic strategy detected: {has_fluorinated_aromatic}")
    return has_fluorinated_aromatic
