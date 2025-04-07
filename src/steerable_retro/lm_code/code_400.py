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
    This function detects if the synthetic route involves a sequence of functional group
    interconversions, particularly acid→amide→isocyanate.
    """
    # Track functional group transformations with molecule identifiers
    acid_to_amide_transformations = []
    amide_to_isocyanate_transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # In retrosynthetic traversal, we're going from product to reactants
            # So for acid→amide→isocyanate sequence, we need to check:
            # 1. amide in product, acid in reactants (amide → acid in retro)
            # 2. isocyanate in product, amide in reactants (isocyanate → amide in retro)

            # Check for amide in product (acid to amide transformation in forward direction)
            if checker.check_fg("Carboxylic acid", product_smiles) and any(
                checker.check_fg("Primary amide", r)
                or checker.check_fg("Secondary amide", r)
                or checker.check_fg("Tertiary amide", r)
                for r in reactants_smiles
            ):
                print(f"Depth {depth}: Acid to amide transformation detected")
                acid_to_amide_transformations.append((product_smiles, depth))

            # Check for isocyanate in product (amide to isocyanate in forward direction)
            if checker.check_fg("Isocyanate", product_smiles) and any(
                checker.check_fg("Primary amide", r)
                or checker.check_fg("Secondary amide", r)
                or checker.check_fg("Tertiary amide", r)
                for r in reactants_smiles
            ):
                print(f"Depth {depth}: Amide to isocyanate transformation detected")
                amide_to_isocyanate_transformations.append((product_smiles, depth))

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if we have the acid→amide→isocyanate sequence
    # For a valid sequence, we need both transformations and they must occur in the correct order
    has_sequence = False

    if acid_to_amide_transformations and amide_to_isocyanate_transformations:
        print(
            f"Found {len(acid_to_amide_transformations)} acid→amide and {len(amide_to_isocyanate_transformations)} amide→isocyanate transformations"
        )

        # Sort transformations by depth to check sequence
        acid_to_amide_transformations.sort(key=lambda x: x[1])
        amide_to_isocyanate_transformations.sort(key=lambda x: x[1])

        # Check if the sequence occurs in the correct order
        for acid_amide_mol, acid_amide_depth in acid_to_amide_transformations:
            for iso_mol, iso_depth in amide_to_isocyanate_transformations:
                # The acid→amide transformation should occur before amide→isocyanate
                # In retrosynthetic traversal, earlier steps have higher depth
                if acid_amide_depth > iso_depth:
                    print(
                        f"Found sequence: acid→amide at depth {acid_amide_depth}, amide→isocyanate at depth {iso_depth}"
                    )
                    has_sequence = True
                    break
            if has_sequence:
                break

    if has_sequence:
        print("Functional group interconversion sequence detected: acid→amide→isocyanate")
    else:
        print("No functional group interconversion sequence detected")

    return has_sequence
