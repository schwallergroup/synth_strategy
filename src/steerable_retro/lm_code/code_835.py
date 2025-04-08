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
    This function detects if a chlorinated pyridine moiety is preserved throughout the synthesis.
    """
    # Track molecules in the main synthetic pathway
    main_pathway_mols = []

    def collect_main_pathway_mols(node, depth=0):
        """Collect molecules that are part of the main synthetic pathway"""
        if node["type"] == "mol":
            # Skip starting materials (in_stock)
            if not node.get("in_stock", False):
                main_pathway_mols.append((node["smiles"], depth))

        # Process children
        for child in node.get("children", []):
            collect_main_pathway_mols(child, depth + 1)

    # Collect main pathway molecules
    collect_main_pathway_mols(route)

    # If no molecules in main pathway, return False
    if not main_pathway_mols:
        print("No molecules found in main synthetic pathway")
        return False

    # Sort by depth (to check from final product to intermediates)
    main_pathway_mols.sort(key=lambda x: x[1])

    # Check if the final product has a chlorinated pyridine
    final_product_smiles = main_pathway_mols[0][0]
    has_pyridine = checker.check_ring("pyridine", final_product_smiles)
    has_aromatic_halide = checker.check_fg("Aromatic halide", final_product_smiles)

    if not (has_pyridine and has_aromatic_halide):
        print(f"Final product does not have chlorinated pyridine: {final_product_smiles}")
        return False

    # Check if all molecules in the main pathway have a chlorinated pyridine
    for mol_smiles, depth in main_pathway_mols:
        has_pyridine = checker.check_ring("pyridine", mol_smiles)
        has_aromatic_halide = checker.check_fg("Aromatic halide", mol_smiles)

        if has_pyridine and has_aromatic_halide:
            print(f"Chlorinated pyridine found in: {mol_smiles}")
        else:
            print(f"Molecule without chlorinated pyridine: {mol_smiles}")
            return False

    print("Chlorinated pyridine preserved throughout the main synthetic pathway")
    return True
