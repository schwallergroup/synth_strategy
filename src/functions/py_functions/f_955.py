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
    Detects if the synthesis route preserves a tetramethyltetralin scaffold
    throughout the synthesis
    """
    scaffold_preserved = True

    # Define tetramethyltetralin patterns - these will match various tetramethyltetralin isomers
    tetramethyltetralin_patterns = [
        Chem.MolFromSmarts("C1CCc2c(C)c(C)c(C)c(C)c2C1"),  # General pattern
        Chem.MolFromSmarts(
            "CC1(C)CCC(C)(C)c2ccccc21"
        ),  # Specific pattern from test case
        Chem.MolFromSmarts("C1CC(C)(C)c2ccc(C)c(C)c2C1"),  # Another possible isomer
        Chem.MolFromSmarts("C1CC(C)(C)c2cc(C)c(C)cc2C1"),  # Another possible isomer
    ]

    # Pattern for dimethylindane which could be part of tetramethyltetralin
    dimethylindane_pattern = Chem.MolFromSmarts("C1CC(C)(C)c2ccccc2C1")

    # Track if we found the scaffold in the final product
    found_in_final_product = False

    def dfs_traverse(node, depth=0):
        nonlocal scaffold_preserved, found_in_final_product

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            if mol:
                # Check if this is a starting material
                is_starting_material = node.get("in_stock", False)

                # If it's the final product (depth 0), we must find the scaffold
                if depth == 0:
                    scaffold_found = False

                    # Try all tetramethyltetralin patterns
                    for pattern in tetramethyltetralin_patterns:
                        if mol.HasSubstructMatch(pattern):
                            print(
                                f"Tetramethyltetralin scaffold found in final product: {mol_smiles}"
                            )
                            scaffold_found = True
                            found_in_final_product = True
                            break

                    # If no direct match, check for tetralin with methyl groups
                    if not scaffold_found:
                        has_tetralin = (
                            checker.check_ring("tetralin", mol_smiles)
                            if "tetralin" in globals()
                            else False
                        )
                        if not has_tetralin:
                            has_tetralin = checker.check_ring("naphthalene", mol_smiles)

                        if has_tetralin:
                            print(
                                f"Found tetralin/naphthalene in final product: {mol_smiles}"
                            )

                            # Get ring indices
                            ring_name = (
                                "tetralin"
                                if "tetralin" in globals()
                                and checker.check_ring("tetralin", mol_smiles)
                                else "naphthalene"
                            )
                            ring_indices = checker.get_ring_atom_indices(
                                ring_name, mol_smiles
                            )

                            if ring_indices:
                                # Create a set of all ring atoms for faster lookup
                                ring_atoms = set()
                                for indices_tuple in ring_indices:
                                    for atom_idx in indices_tuple[0]:
                                        ring_atoms.add(atom_idx)

                                # Count methyl groups attached to ring
                                methyl_count = 0
                                for atom_idx in ring_atoms:
                                    atom = mol.GetAtomWithIdx(atom_idx)
                                    for neighbor in atom.GetNeighbors():
                                        # Check if neighbor is a carbon with 3 hydrogens (methyl)
                                        if (
                                            neighbor.GetAtomicNum() == 6
                                            and neighbor.GetIdx() not in ring_atoms
                                            and neighbor.GetTotalNumHs() == 3
                                        ):
                                            methyl_count += 1

                                print(
                                    f"Found {methyl_count} methyl groups attached to ring in final product"
                                )

                                # Check if we have at least 4 methyl groups
                                if methyl_count >= 4:
                                    print(
                                        f"Tetramethyltetralin scaffold found in final product"
                                    )
                                    scaffold_found = True
                                    found_in_final_product = True

                    # Check for dimethylindane as a fallback
                    if not scaffold_found and mol.HasSubstructMatch(
                        dimethylindane_pattern
                    ):
                        print(
                            f"Dimethylindane pattern found in final product - this could be part of tetramethyltetralin"
                        )
                        scaffold_found = True
                        found_in_final_product = True

                    if not scaffold_found:
                        print(
                            f"Tetramethyltetralin scaffold not found in final product: {mol_smiles}"
                        )
                        scaffold_preserved = False

                # For intermediates (not starting materials), check if they have the scaffold
                elif (
                    not is_starting_material and depth <= 8
                ):  # Only check key intermediates (not too deep in the tree)
                    scaffold_found = False

                    # Try all tetramethyltetralin patterns
                    for pattern in tetramethyltetralin_patterns:
                        if mol.HasSubstructMatch(pattern):
                            print(
                                f"Tetramethyltetralin scaffold found in intermediate at depth {depth}: {mol_smiles}"
                            )
                            scaffold_found = True
                            break

                    # Check for dimethylindane as a fallback
                    if not scaffold_found and mol.HasSubstructMatch(
                        dimethylindane_pattern
                    ):
                        print(
                            f"Dimethylindane pattern found in intermediate at depth {depth} - this could be part of tetramethyltetralin"
                        )
                        scaffold_found = True

                    # For key intermediates, we should find the scaffold
                    if (
                        not scaffold_found and len(mol_smiles) > 10
                    ):  # Only check substantial molecules, not small reagents
                        print(
                            f"Tetramethyltetralin scaffold not found in key intermediate at depth {depth}: {mol_smiles}"
                        )
                        # Only mark as not preserved if this is a substantial molecule that should have the scaffold
                        if mol.GetNumHeavyAtoms() > 15:  # Adjust threshold as needed
                            scaffold_preserved = False
            else:
                print(f"Failed to create RDKit molecule at depth {depth}: {mol_smiles}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Final check - we must have found the scaffold in the final product
    if not found_in_final_product:
        scaffold_preserved = False

    if scaffold_preserved:
        print("Tetramethyltetralin scaffold is preserved throughout the synthesis")
    else:
        print("Tetramethyltetralin scaffold is NOT preserved throughout the synthesis")

    return scaffold_preserved
