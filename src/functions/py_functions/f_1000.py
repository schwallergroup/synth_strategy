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
    Detects if the synthesis involves a nucleophilic aromatic substitution on a pyrimidine ring.
    """
    nas_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nas_found

        if node["type"] == "reaction" and depth <= 1:  # Late in synthesis (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check if any reactant has a pyrimidine ring with a halide
                    reactant_has_pyrimidine_halide = False
                    for reactant in reactants:
                        if checker.check_ring(
                            "pyrimidine", reactant
                        ) and checker.check_fg("Aromatic halide", reactant):
                            print(
                                f"Found pyrimidine with aromatic halide in reactant: {reactant}"
                            )
                            reactant_has_pyrimidine_halide = True
                            break

                    # Check if product has a pyrimidine ring
                    product_has_pyrimidine = checker.check_ring("pyrimidine", product)

                    # Check if this is a nucleophilic aromatic substitution reaction
                    is_nas_reaction = (
                        checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                        or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                        or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    )

                    # Additional check: if reaction types don't match, look for amine/nitrogen nucleophile in product
                    if (
                        not is_nas_reaction
                        and reactant_has_pyrimidine_halide
                        and product_has_pyrimidine
                    ):
                        # Check if product has an amine group where the halide was
                        if (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                            or checker.check_fg("Aniline", product)
                        ):
                            print(
                                f"Found potential NAS based on structural changes at depth {depth}"
                            )
                            is_nas_reaction = True

                    if (
                        reactant_has_pyrimidine_halide
                        and product_has_pyrimidine
                        and is_nas_reaction
                    ):
                        print(
                            f"Nucleophilic aromatic substitution on pyrimidine detected at depth {depth}"
                        )
                        nas_found = True
                    else:
                        if not reactant_has_pyrimidine_halide:
                            print(
                                f"No pyrimidine with halide found in reactants at depth {depth}"
                            )
                        if not product_has_pyrimidine:
                            print(f"No pyrimidine found in product at depth {depth}")
                        if not is_nas_reaction:
                            print(
                                f"Reaction at depth {depth} is not a nucleophilic aromatic substitution"
                            )

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nas_found
