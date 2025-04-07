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
    Detects if the synthesis involves a product with multiple rings.
    """
    has_multiple_rings = False
    ring_types_found = set()

    def dfs_traverse(node, depth=0):
        nonlocal has_multiple_rings, ring_types_found

        if node["type"] == "mol":
            # Check if this is a non-stock molecule (potential product)
            if "in_stock" not in node or not node["in_stock"]:
                smiles = node["smiles"]
                mol = Chem.MolFromSmiles(smiles)

                if mol:
                    # Get ring information
                    ring_info = mol.GetRingInfo()
                    num_rings = ring_info.NumRings()

                    # Check for specific ring types
                    ring_types = []
                    for ring_name in [
                        "benzene",
                        "pyridine",
                        "pyrrole",
                        "furan",
                        "thiophene",
                        "pyrimidine",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "piperidine",
                        "tetrahydrofuran",
                        "cyclohexane",
                        "cyclopentane",
                    ]:
                        if checker.check_ring(ring_name, smiles):
                            ring_types.append(ring_name)
                            ring_types_found.add(ring_name)

                    # If we have multiple rings, set the flag
                    if num_rings >= 2:
                        has_multiple_rings = True
                        print(
                            f"Molecule at depth {depth} contains {num_rings} rings: {', '.join(ring_types)}"
                        )

        # Traverse children (in retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_multiple_rings:
        print(
            f"Multi-cyclic structure strategy: True. Ring types found: {', '.join(ring_types_found)}"
        )
    else:
        print(f"Multi-cyclic structure strategy: False. No multi-cyclic structures detected.")

    return has_multiple_rings
