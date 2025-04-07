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
    Detects if the synthesis preserves a 2,4-disubstituted aromatic core
    (with F and Cl substituents) throughout the synthesis.
    """
    # Track if we've found the core in the final product
    core_found_in_final = False
    # Track if the core is preserved throughout the main synthetic pathway
    core_preserved = True

    # First check if the final product has the core
    final_product = route
    if final_product["type"] == "mol":
        mol_smiles = final_product["smiles"]
        has_benzene = checker.check_ring("benzene", mol_smiles)
        has_f = checker.check_fg("Aromatic halide", mol_smiles) and "F" in mol_smiles
        has_cl = checker.check_fg("Aromatic halide", mol_smiles) and "Cl" in mol_smiles

        if has_benzene and has_f and has_cl:
            core_found_in_final = True
            print(f"Core found in final product: {mol_smiles}")
        else:
            print(f"Core not found in final product: {mol_smiles}")
            return False  # If final product doesn't have core, return False immediately

    # Track the main synthetic pathway
    main_pathway_nodes = []

    def trace_main_pathway(node, depth=0):
        """Trace the main synthetic pathway by following product to reactants"""
        if node["type"] == "mol":
            main_pathway_nodes.append((node, depth))

            # For molecule nodes, check all reaction children
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    trace_main_pathway(child, depth + 1)

        elif node["type"] == "reaction":
            # For reaction nodes, only follow the main product path
            # (the first molecule child, which is typically the main reactant)
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    trace_main_pathway(child, depth + 1)
                    break  # Only follow the first non-stock molecule

    # Start tracing from the final product
    trace_main_pathway(route)

    # Sort by depth (later stages first)
    main_pathway_nodes.sort(key=lambda x: x[1])

    # Check each molecule in the main pathway
    for node, depth in main_pathway_nodes:
        mol_smiles = node["smiles"]
        has_benzene = checker.check_ring("benzene", mol_smiles)
        has_f = checker.check_fg("Aromatic halide", mol_smiles) and "F" in mol_smiles
        has_cl = checker.check_fg("Aromatic halide", mol_smiles) and "Cl" in mol_smiles

        if has_benzene and has_f and has_cl:
            print(f"Core found in molecule at depth {depth}: {mol_smiles}")
        else:
            # If we're past the first molecule and the core is missing, mark as not preserved
            if depth > 0:
                print(f"Core not found in molecule at depth {depth}: {mol_smiles}")
                core_preserved = False

    return core_preserved and core_found_in_final
