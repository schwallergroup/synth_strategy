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
    This function detects if the synthesis involves connecting a fluorinated aromatic
    compound with a sulfonamide containing a cyclic amine.
    """
    has_fluorinated_aromatic = False
    has_cyclic_amine = False
    has_connected_structure = False

    # List of cyclic amines to check
    cyclic_amine_rings = [
        "aziridine",
        "azetidine",
        "pyrrolidine",
        "piperidine",
        "azepane",
        "diazepane",
        "morpholine",
        "piperazine",
        "thiomorpholine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_fluorinated_aromatic, has_cyclic_amine, has_connected_structure

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for fluorinated aromatic
            if checker.check_fg("Aromatic halide", mol_smiles) and "F" in mol_smiles:
                has_fluorinated_aromatic = True
                print(f"Detected fluorinated aromatic compound: {mol_smiles}")

            # Check for cyclic amine
            for ring in cyclic_amine_rings:
                if checker.check_ring(ring, mol_smiles):
                    has_cyclic_amine = True
                    print(f"Detected cyclic amine ({ring}): {mol_smiles}")
                    break

            # Check if this molecule already has both components connected via sulfonamide
            if (
                checker.check_fg("Aromatic halide", mol_smiles)
                and "F" in mol_smiles
                and checker.check_fg("Sulfonamide", mol_smiles)
                and any(checker.check_ring(ring, mol_smiles) for ring in cyclic_amine_rings)
            ):
                has_connected_structure = True
                print(
                    f"Detected molecule with fluorinated aromatic and cyclic amine connected via sulfonamide: {mol_smiles}"
                )

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains both fluorinated aromatic, sulfonamide, and cyclic amine
                product_has_fluoro = checker.check_fg("Aromatic halide", product) and "F" in product
                product_has_sulfonamide = checker.check_fg("Sulfonamide", product)
                product_has_cyclic_amine = any(
                    checker.check_ring(ring, product) for ring in cyclic_amine_rings
                )

                if product_has_fluoro and product_has_sulfonamide and product_has_cyclic_amine:
                    # Check if the product has all three components
                    has_connected_structure = True
                    print(
                        f"Detected product with fluorinated aromatic and cyclic amine connected via sulfonamide: {product}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if all conditions are met
    result = has_fluorinated_aromatic and has_cyclic_amine and has_connected_structure
    print(f"Final result: {result}")
    print(f"- Has fluorinated aromatic: {has_fluorinated_aromatic}")
    print(f"- Has cyclic amine: {has_cyclic_amine}")
    print(f"- Has connected structure: {has_connected_structure}")

    return result
