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
    Detects if the synthetic route involves preservation of fluorinated aromatic groups
    throughout the synthesis.
    """
    # Track fluorinated molecules at each depth
    fluorinated_at_depth = {}

    def dfs_traverse(node, depth=0):
        # Process molecule nodes
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for fluorinated aromatics (either aryl-F or CF3)
            has_aryl_f = checker.check_fg(
                "Aromatic halide", mol_smiles
            ) or checker.check_fg("Trifluoro group", mol_smiles)

            if has_aryl_f:
                print(
                    f"Detected fluorinated aryl at depth {depth}, SMILES: {mol_smiles}"
                )
                if depth not in fluorinated_at_depth:
                    fluorinated_at_depth[depth] = []
                fluorinated_at_depth[depth].append(mol_smiles)

        # Process reaction nodes
        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # Check if fluorinated groups are preserved in this reaction
            product_has_fluorine = checker.check_fg(
                "Aromatic halide", product
            ) or checker.check_fg("Trifluoro group", product)
            reactants_have_fluorine = any(
                checker.check_fg("Aromatic halide", r)
                or checker.check_fg("Trifluoro group", r)
                for r in reactants
            )

            if product_has_fluorine and reactants_have_fluorine:
                print(
                    f"Fluorinated group preserved in reaction at depth {depth}, RSMI: {rsmi}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have fluorinated groups at multiple depths
    if len(fluorinated_at_depth) >= 2:
        depths = sorted(fluorinated_at_depth.keys())
        depth_range = max(depths) - min(depths)

        # Check if the depth range is significant (at least 2 levels)
        if depth_range >= 2:
            print(
                f"Found fluorinated groups preserved across multiple depths: {min(depths)}-{max(depths)}"
            )

            # Check if we have fluorinated groups at both early and late stages
            early_stage = max(depths)
            late_stage = min(depths)

            if early_stage - late_stage >= 2:
                print(
                    f"Fluorinated groups preserved from early stage (depth {early_stage}) to late stage (depth {late_stage})"
                )
                return True

    return False
