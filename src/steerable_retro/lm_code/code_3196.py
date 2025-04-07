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
    This function detects a strategy involving sequential SNAr reactions
    on a pyrazole scaffold with late-stage ether formation.
    """
    # Track key features
    has_pyrazole = False
    snAr_reactions = 0
    late_stage_ether_formation = False
    has_halogen_intermediates = False

    def dfs_traverse(node, depth=0):
        nonlocal has_pyrazole, snAr_reactions, late_stage_ether_formation, has_halogen_intermediates

        if node["type"] == "mol":
            # Check for pyrazole scaffold
            if node["smiles"]:
                # Check for pyrazole ring using both the checker and a pattern-based approach
                mol_smiles = node["smiles"]
                if checker.check_ring("pyrazole", mol_smiles):
                    has_pyrazole = True
                    print(f"Pyrazole scaffold detected via checker in: {mol_smiles}")
                # Backup check for pyrazole pattern in SMILES
                elif "c1cnc" in mol_smiles and "n1" in mol_smiles:
                    has_pyrazole = True
                    print(f"Pyrazole scaffold detected via pattern in: {mol_smiles}")
                # Another common pyrazole pattern
                elif "c1nc" in mol_smiles and "cn1" in mol_smiles:
                    has_pyrazole = True
                    print(f"Pyrazole scaffold detected via alternate pattern in: {mol_smiles}")

                # Check for halogen intermediates using functional group checker
                if checker.check_fg("Aromatic halide", mol_smiles):
                    has_halogen_intermediates = True
                    print(f"Aromatic halide detected in molecule: {mol_smiles}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr reactions using reaction checkers
                is_snar = (
                    checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                )

                # If not directly identified, check for characteristic patterns
                if not is_snar:
                    # Check if any reactant has aromatic halide
                    has_aromatic_halide = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants if r
                    )

                    # Check if product has new ether or amine where halogen was
                    has_new_bond = (
                        checker.check_fg("Ether", product)
                        or checker.check_fg("Aniline", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    )

                    is_snar = has_aromatic_halide and has_new_bond

                if is_snar:
                    snAr_reactions += 1
                    print(f"SNAr reaction detected at depth {depth}: {rsmi}")

                    # Check for late-stage ether formation (depth 0 or 1)
                    if depth <= 1 and checker.check_fg("Ether", product):
                        late_stage_ether_formation = True
                        print(f"Late-stage ether formation detected at depth {depth}: {product}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if this strategy is present
    # Note: has_halogen_intermediates is redundant with SNAr reactions requirement
    strategy_present = has_pyrazole and snAr_reactions >= 2 and late_stage_ether_formation

    print(f"Strategy detection results:")
    print(f"  Pyrazole scaffold: {has_pyrazole}")
    print(f"  SNAr reactions: {snAr_reactions}")
    print(f"  Late-stage ether formation: {late_stage_ether_formation}")
    print(f"  Halogen intermediates: {has_halogen_intermediates}")
    print(f"  Strategy present: {strategy_present}")

    return strategy_present
