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
    This function detects if the synthesis incorporates multiple heteroaromatic rings.
    """
    # List of common heteroaromatic rings to check
    heteroaromatic_rings = [
        "pyridine",
        "pyrazole",
        "imidazole",
        "thiazole",
        "oxazole",
        "furan",
        "thiophene",
        "pyrimidine",
        "pyrazine",
        "indole",
        "benzimidazole",
        "benzothiazole",
        "benzoxazole",
        "quinoline",
        "isoquinoline",
        "triazole",
        "tetrazole",
        "isoxazole",
        "isothiazole",
    ]

    # Use a set to track found heteroaromatic rings
    found_rings = set()

    def dfs_traverse(node):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for each heteroaromatic ring type
            for ring_type in heteroaromatic_rings:
                if (
                    checker.check_ring(ring_type, mol_smiles)
                    and ring_type not in found_rings
                ):
                    found_rings.add(ring_type)
                    print(f"Found {ring_type} ring in molecule: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Found heteroaromatic rings: {found_rings}")
    print(f"Total unique heteroaromatic rings found: {len(found_rings)}")

    # Return True if at least 2 different heteroaromatic rings are found
    return len(found_rings) >= 2
