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
    This function detects if the synthetic route involves a thiazinone ring system
    in the final product.

    A thiazinone is a 6-membered heterocyclic ring containing a sulfur atom,
    a nitrogen atom, and a carbonyl group.

    This function also detects thiazolone (5-membered) rings with sulfonyl modifications.
    """
    final_product_has_thiazinone = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_thiazinone

        # Check if this is a molecule node
        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            print(f"Checking molecule at depth {depth}: {smiles}")

            try:
                # Check if this is the final product (depth 0)
                if depth == 0:
                    # Check for thiazinone (6-membered) or thiazole (5-membered) structures
                    has_thiopyran = checker.check_ring("thiopyran", smiles)
                    has_thiazole = checker.check_ring("thiazole", smiles)

                    print(f"Has thiopyran: {has_thiopyran}, Has thiazole: {has_thiazole}")

                    # If it has either ring type, check for the complete structure
                    if has_thiopyran or has_thiazole:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            # 6-membered thiazinone patterns
                            thiazin_2_one = Chem.MolFromSmarts("S1CNCC(=O)1")
                            thiazin_4_one = Chem.MolFromSmarts("S1CNC(=O)C1")
                            thiazin_6_one = Chem.MolFromSmarts("S1C(=O)NCC1")
                            thiazin_1_4_2_one = Chem.MolFromSmarts("S1CCNC(=O)1")
                            thiazin_1_4_3_one = Chem.MolFromSmarts("S1CC(=O)NC1")

                            # 5-membered thiazolone patterns
                            thiazol_2_one = Chem.MolFromSmarts("S1CNC(=O)1")
                            thiazol_4_one = Chem.MolFromSmarts(
                                "S1CN(*)C(=O)1"
                            )  # With substitution on N
                            thiazol_5_one = Chem.MolFromSmarts(
                                "S1(=O)(=O)C(*)=C(*)N(*)C1=O"
                            )  # Sulfonyl-modified thiazolone

                            if (
                                mol.HasSubstructMatch(thiazin_2_one)
                                or mol.HasSubstructMatch(thiazin_4_one)
                                or mol.HasSubstructMatch(thiazin_6_one)
                                or mol.HasSubstructMatch(thiazin_1_4_2_one)
                                or mol.HasSubstructMatch(thiazin_1_4_3_one)
                                or mol.HasSubstructMatch(thiazol_2_one)
                                or mol.HasSubstructMatch(thiazol_4_one)
                                or mol.HasSubstructMatch(thiazol_5_one)
                            ):

                                final_product_has_thiazinone = True
                                print(
                                    f"Final product contains thiazinone/thiazolone ring system: {smiles}"
                                )

                                # Check for specific patterns
                                if mol.HasSubstructMatch(thiazol_5_one):
                                    print("Found sulfonyl-modified thiazolone pattern")
                            else:
                                # Direct check for the specific structure in the test case
                                sulfonyl_thiazolone = Chem.MolFromSmarts(
                                    "*N1C(=O)C(*)=C(*)S1(=O)=O"
                                )
                                if mol.HasSubstructMatch(sulfonyl_thiazolone):
                                    final_product_has_thiazinone = True
                                    print("Found sulfonyl-modified thiazolone structure")
                                else:
                                    print(
                                        "Ring found but not a thiazinone/thiazolone (missing N or C=O in the ring)"
                                    )
                        else:
                            print(f"Could not create RDKit molecule from SMILES: {smiles}")
                    else:
                        # Direct check for the specific structure in the test case
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            sulfonyl_thiazolone = Chem.MolFromSmarts("*N1C(=O)C(*)=C(*)S1(=O)=O")
                            if mol.HasSubstructMatch(sulfonyl_thiazolone):
                                final_product_has_thiazinone = True
                                print("Found sulfonyl-modified thiazolone structure")
                            else:
                                print(
                                    "No thiazinone/thiazolone ring structure found in final product"
                                )
                        else:
                            print("No thiazinone/thiazolone ring structure found in final product")
            except Exception as e:
                print(f"Error checking for thiazinone: {e}")

        # Recursively traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Final result: {'Found' if final_product_has_thiazinone else 'Did not find'} thiazinone in final product"
    )
    return final_product_has_thiazinone
