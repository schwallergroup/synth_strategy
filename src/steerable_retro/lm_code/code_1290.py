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
    Detects a strategy where a benzofuran core scaffold is maintained throughout the synthesis
    while peripheral functional groups are modified.
    """
    # Track if benzofuran is present in all molecules
    all_have_benzofuran = True
    molecule_count = 0

    def dfs_traverse(node):
        nonlocal all_have_benzofuran, molecule_count

        if node["type"] == "mol" and not node.get("in_stock", False):
            molecule_count += 1
            smiles = node["smiles"]

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    print(f"Failed to parse molecule SMILES: {smiles}")
                    all_have_benzofuran = False
                    return

                # Check for benzofuran core using the checker function
                has_benzofuran = checker.check_ring("benzofuran", smiles)

                # Additional check with a more general pattern for benzofuran
                if not has_benzofuran:
                    # Try a more general pattern for benzofuran (benzene fused to furan)
                    benzofuran_pattern = Chem.MolFromSmarts("c1ccc2c(c1)cco2")
                    if mol.HasSubstructMatch(benzofuran_pattern):
                        has_benzofuran = True
                        print(f"SMARTS pattern detected benzofuran in: {smiles}")

                if not has_benzofuran:
                    all_have_benzofuran = False
                    print(f"Molecule without benzofuran core: {smiles}")
                else:
                    print(f"Molecule with benzofuran core: {smiles}")
            except Exception as e:
                print(f"Error processing molecule SMILES: {smiles} - {e}")
                all_have_benzofuran = False

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if all non-starting molecules contain benzofuran and we have multiple molecules
    print(f"All have benzofuran: {all_have_benzofuran}, Molecule count: {molecule_count}")
    return all_have_benzofuran and molecule_count > 1
