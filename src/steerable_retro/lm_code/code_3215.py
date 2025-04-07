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
    Detects if the synthesis route involves heterocyclic scaffolds,
    specifically piperazine and tetrahydroisoquinoline motifs.
    """
    from rdkit import Chem

    has_piperazine = False
    has_tetrahydroisoquinoline = False

    # Define SMARTS patterns for more precise matching
    piperazinone_pattern = Chem.MolFromSmarts("N1CCN(*)CC(=O)1")

    def dfs_traverse(node, depth=0):
        nonlocal has_piperazine, has_tetrahydroisoquinoline

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            print(f"Checking molecule at depth {depth}: {mol_smiles}")

            # Check for piperazine using the checker function
            if checker.check_ring("piperazine", mol_smiles):
                print(f"Found piperazine scaffold at depth {depth} in molecule: {mol_smiles}")
                has_piperazine = True

            # Check for piperazine-like structures (including piperazinone)
            elif "N1CCN" in mol_smiles:
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol and mol.HasSubstructMatch(piperazinone_pattern):
                    print(f"Found piperazinone scaffold at depth {depth} in molecule: {mol_smiles}")
                    has_piperazine = True

            # Check for tetrahydroisoquinoline
            # Direct check for tetrahydroisoquinoline structure
            if checker.check_ring("isoquinoline", mol_smiles):
                # Verify it's the tetrahydro form
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Check if it contains a 1,2,3,4-tetrahydroisoquinoline structure
                    # This is a simplified check - in a real scenario, we would need more sophisticated analysis
                    if "CCc2ccc" in mol_smiles and "cc2C1" in mol_smiles:
                        print(
                            f"Found tetrahydroisoquinoline-like scaffold at depth {depth} in molecule: {mol_smiles}"
                        )
                        has_tetrahydroisoquinoline = True

            # Additional check for tetrahydroisoquinoline-like structures
            elif "CCc2ccc" in mol_smiles and "cc2C1" in mol_smiles:
                print(
                    f"Found potential tetrahydroisoquinoline scaffold at depth {depth} in molecule: {mol_smiles}"
                )
                has_tetrahydroisoquinoline = True

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rxn_smiles}")

            # Extract product from reaction SMILES
            product = rxn_smiles.split(">")[-1]

            # Check if the product contains our target scaffolds
            if checker.check_ring("piperazine", product):
                print(f"Found piperazine scaffold in reaction product at depth {depth}")
                has_piperazine = True

            # Check for piperazine-like structures in product
            elif "N1CCN" in product:
                mol = Chem.MolFromSmiles(product)
                if mol and mol.HasSubstructMatch(piperazinone_pattern):
                    print(f"Found piperazinone scaffold in reaction product at depth {depth}")
                    has_piperazine = True

            # Check for tetrahydroisoquinoline in product
            if checker.check_ring("isoquinoline", product):
                mol = Chem.MolFromSmiles(product)
                if mol:
                    # Check if it contains a 1,2,3,4-tetrahydroisoquinoline structure
                    if "CCc2ccc" in product and "cc2C1" in product:
                        print(
                            f"Found tetrahydroisoquinoline-like scaffold in reaction product at depth {depth}"
                        )
                        has_tetrahydroisoquinoline = True

            # Additional check for tetrahydroisoquinoline-like structures in product
            elif "CCc2ccc" in product and "cc2C1" in product:
                print(
                    f"Found potential tetrahydroisoquinoline scaffold in reaction product at depth {depth}"
                )
                has_tetrahydroisoquinoline = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if either scaffold is found
    result = has_piperazine or has_tetrahydroisoquinoline
    print(
        f"Final result: {result} (Piperazine: {has_piperazine}, Tetrahydroisoquinoline: {has_tetrahydroisoquinoline})"
    )
    return result
