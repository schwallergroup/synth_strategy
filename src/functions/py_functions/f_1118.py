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
    Detects if a fluorinated aromatic group is introduced in a late stage of the synthesis.
    Late stage is defined as the first half of the synthesis (lower depth in retrosynthesis).
    """
    max_depth = 0
    fluorinated_introduction_depth = None

    def has_fluorinated_aromatic(smiles):
        """Check if molecule contains a fluorinated aromatic group"""
        # Check for aromatic halide (which includes fluorine) or trifluoro group
        return checker.check_fg("Aromatic halide", smiles) or checker.check_fg(
            "Trifluoro group", smiles
        )

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, fluorinated_introduction_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            try:
                if "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0]
                    product = rsmi.split(">")[2]  # Get the product part

                    print(f"Analyzing reaction at depth {depth}")
                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check if fluorinated group is in product but not in all reactants
                    product_has_fluorine = has_fluorinated_aromatic(product)
                    reactant_mols = reactants.split(".")
                    all_reactants_have_fluorine = all(
                        has_fluorinated_aromatic(r) for r in reactant_mols
                    )

                    print(f"Product has fluorinated group: {product_has_fluorine}")
                    print(
                        f"All reactants have fluorinated group: {all_reactants_have_fluorine}"
                    )

                    if product_has_fluorine and not all_reactants_have_fluorine:
                        # This reaction introduces a fluorinated group
                        if (
                            fluorinated_introduction_depth is None
                            or depth < fluorinated_introduction_depth
                        ):
                            fluorinated_introduction_depth = depth
                            print(
                                f"Found fluorinated group introduction at depth {depth}"
                            )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if fluorinated group is introduced in the first half of the synthesis
    if fluorinated_introduction_depth is not None:
        is_late_stage = fluorinated_introduction_depth <= max_depth / 2
        print(f"Max depth: {max_depth}")
        print(f"Fluorinated introduction depth: {fluorinated_introduction_depth}")
        print(f"Fluorinated group introduction is late stage: {is_late_stage}")
        return is_late_stage
    else:
        print("No fluorinated group introduction found")

    return False
