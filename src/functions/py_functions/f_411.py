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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects a strategy involving THP protection of an alcohol
    followed by deprotection later in the synthesis.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_found, deprotection_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthetic traversal:
                # - For protection: product has alcohol, reactants have THP-protected alcohol
                # - For deprotection: product has THP-protected alcohol, reactants have alcohol

                # Check for THP protection (in retrosynthesis: alcohol ← protected alcohol)
                if not protection_found:
                    # Check if product has alcohol
                    product_has_alcohol = any(
                        checker.check_fg("Primary alcohol", p)
                        or checker.check_fg("Secondary alcohol", p)
                        or checker.check_fg("Tertiary alcohol", p)
                        for p in product.split(".")
                        if p
                    )

                    # Check if reactants have THP-protected structure
                    reactants_have_thp_protected = False
                    for r in reactants:
                        if r and checker.check_ring("tetrahydropyran", r):
                            # Check if this THP is connected to an oxygen (ether)
                            mol = Chem.MolFromSmiles(r)
                            if mol:
                                for atom in mol.GetAtoms():
                                    if (
                                        atom.GetSymbol() == "O"
                                        and atom.GetDegree() == 2
                                    ):
                                        # Check if one neighbor is in THP ring and other is carbon
                                        neighbors = [
                                            mol.GetAtomWithIdx(n.GetIdx())
                                            for n in atom.GetNeighbors()
                                        ]
                                        if any(n.IsInRing() for n in neighbors) and any(
                                            n.GetSymbol() == "C" and not n.IsInRing()
                                            for n in neighbors
                                        ):
                                            reactants_have_thp_protected = True
                                            break

                    if product_has_alcohol and reactants_have_thp_protected:
                        protection_found = True
                        print(f"THP protection detected at depth {depth}")

                # Check for THP deprotection (in retrosynthesis: protected alcohol ← alcohol)
                if not deprotection_found:
                    # Check if product has THP-protected structure
                    product_has_thp_protected = False
                    if checker.check_ring("tetrahydropyran", product):
                        mol = Chem.MolFromSmiles(product)
                        if mol:
                            for atom in mol.GetAtoms():
                                if atom.GetSymbol() == "O" and atom.GetDegree() == 2:
                                    # Check if one neighbor is in THP ring and other is carbon
                                    neighbors = [
                                        mol.GetAtomWithIdx(n.GetIdx())
                                        for n in atom.GetNeighbors()
                                    ]
                                    if any(n.IsInRing() for n in neighbors) and any(
                                        n.GetSymbol() == "C" and not n.IsInRing()
                                        for n in neighbors
                                    ):
                                        product_has_thp_protected = True
                                        break

                    # Check if reactants have alcohol
                    reactants_have_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                        if r
                    )

                    if product_has_thp_protected and reactants_have_alcohol:
                        deprotection_found = True
                        print(f"THP deprotection detected at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return protection_found and deprotection_found
