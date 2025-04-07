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
    This function detects aromatic fragment coupling via nucleophilic aromatic substitution.
    """
    found_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]

                # Check for nucleophilic aromatic substitution reactions directly
                if (
                    checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                ):
                    found_snar = True
                    print(
                        f"Nucleophilic aromatic substitution detected at depth {depth} with reaction type check"
                    )
                    return

                # If direct reaction check fails, check for characteristic patterns
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for aromatic halide in reactants
                has_aromatic_halide = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                )

                if has_aromatic_halide:
                    # Check for nitrogen nucleophiles in reactants
                    has_nitrogen_nucleophile = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        or checker.check_fg("Azide", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants_smiles
                    )

                    if has_nitrogen_nucleophile:
                        # Check if product has new N-aromatic connection
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol:
                            # If we have both aromatic halide and nitrogen nucleophile in reactants,
                            # and the product is formed, it's likely a nucleophilic aromatic substitution
                            found_snar = True
                            print(
                                f"Nucleophilic aromatic substitution detected at depth {depth} with FG checks"
                            )
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_snar
