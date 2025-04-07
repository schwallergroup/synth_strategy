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
    Detects if the synthesis maintains a complex bicyclic core with chlorinated aromatic rings.
    """
    steps_with_core = 0
    total_steps = 0

    def has_bicyclic_core(mol_smiles):
        """Check if molecule has the required bicyclic core with chlorinated aromatic rings"""
        # Check for benzene rings
        has_benzene = checker.check_ring("benzene", mol_smiles)

        # Check for chlorinated aromatic
        has_aromatic_halide = checker.check_fg("Aromatic halide", mol_smiles)

        # Check for lactam structure (contains both cyclic amide and benzene)
        has_amide = (
            checker.check_fg("Primary amide", mol_smiles)
            or checker.check_fg("Secondary amide", mol_smiles)
            or checker.check_fg("Tertiary amide", mol_smiles)
        )

        # Check for bicyclic structure
        bicyclic_structure = False
        for ring1 in ["cyclohexane", "cyclopentane", "piperidine", "pyrrolidine"]:
            if checker.check_ring(ring1, mol_smiles):
                bicyclic_structure = True
                break

        return has_benzene and has_aromatic_halide and has_amide and bicyclic_structure

    def check_core_preservation_in_reaction(reaction_node):
        """Check if a reaction preserves the bicyclic core structure"""
        try:
            rsmi = reaction_node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has the core
            product_has_core = has_bicyclic_core(product)

            # If product doesn't have core, no preservation
            if not product_has_core:
                return False

            # Check if at least one reactant has the core
            reactant_has_core = any(has_bicyclic_core(r) for r in reactants)

            # Core preservation means it existed in reactants and continues in product
            return reactant_has_core
        except Exception as e:
            print(f"Error checking core preservation in reaction: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_core, total_steps

        if node["type"] == "reaction":
            # For reaction nodes, check if the core is preserved
            total_steps += 1
            if check_core_preservation_in_reaction(node):
                steps_with_core += 1
                print(f"Core preserved in reaction at depth {depth}")
            else:
                print(f"Core NOT preserved in reaction at depth {depth}")

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Bicyclic core preserved in {steps_with_core}/{total_steps} steps")

    # If the core is preserved in most steps (>80%) and we have at least one step
    if total_steps > 0 and steps_with_core / total_steps > 0.8:
        return True

    return False
