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
    Detects preservation of difluoromethoxy group throughout the synthesis
    """
    # Track if difluoromethoxy group is present in the main synthetic pathway
    all_steps_have_difluoromethoxy = True
    step_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal all_steps_have_difluoromethoxy, step_count

        if node["type"] == "mol" and "smiles" in node:
            # Only check molecules that are part of the main synthetic pathway
            # (not in_stock reagents or catalysts)
            if depth == 0 or (node.get("children") and len(node.get("children")) > 0):
                mol_smiles = node["smiles"]

                # Check for difluoromethoxy group using the checker function
                has_difluoromethoxy = checker.check_fg("Trifluoro group", mol_smiles)

                if not has_difluoromethoxy:
                    # Try alternative check for difluoromethoxy (OCF2H)
                    difluoromethoxy_pattern = Chem.MolFromSmarts("OC(F)F")
                    mol = Chem.MolFromSmiles(mol_smiles)
                    has_difluoromethoxy = mol and mol.HasSubstructMatch(difluoromethoxy_pattern)

                if not has_difluoromethoxy:
                    all_steps_have_difluoromethoxy = False
                    print(f"No difluoromethoxy group found in main pathway molecule: {mol_smiles}")
                else:
                    print(f"Found difluoromethoxy group in main pathway molecule: {mol_smiles}")

                step_count += 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Only return True if we've checked at least one molecule and all have difluoromethoxy
    return all_steps_have_difluoromethoxy and step_count > 0
