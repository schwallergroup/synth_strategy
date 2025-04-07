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
    This function detects if the synthesis maintains a pyridazinone core throughout.
    """
    # Track if all non-starting material molecules have pyridazinone core
    all_have_core = True

    def dfs_traverse(node, depth=0):
        nonlocal all_have_core

        if node["type"] == "mol" and "smiles" in node:
            # Skip checking starting materials
            if node.get("in_stock", False):
                return

            # Check for pyridazinone core (pyridazine ring with carbonyl)
            has_pyridazine = checker.check_ring("pyridazine", node["smiles"])

            if has_pyridazine:
                # Check for carbonyl group attached to pyridazine
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Pyridazinone pattern: pyridazine with carbonyl
                    pyridazinone_pattern = Chem.MolFromSmarts("n1ncc(=O)cc1")
                    alt_pattern = Chem.MolFromSmarts("n1nc(=O)ccc1")

                    if not (
                        mol.HasSubstructMatch(pyridazinone_pattern)
                        or mol.HasSubstructMatch(alt_pattern)
                    ):
                        print(
                            f"Molecule at depth {depth} has pyridazine but no carbonyl: {node['smiles']}"
                        )
                        all_have_core = False
            else:
                print(f"Molecule at depth {depth} missing pyridazine ring: {node['smiles']}")
                all_have_core = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return all_have_core
