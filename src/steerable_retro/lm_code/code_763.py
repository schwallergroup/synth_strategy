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
    This function detects if the synthetic route involves nucleophilic aromatic substitution.
    """
    has_snar = False

    def dfs_traverse(node):
        nonlocal has_snar

        if node["type"] == "reaction" and not has_snar:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for nucleophilic aromatic substitution reactions directly
                if (
                    checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                ):
                    print(f"Found nucleophilic aromatic substitution reaction: {rsmi}")
                    has_snar = True
                else:
                    # If direct reaction check fails, check for characteristic patterns
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for aryl halide in reactants
                    has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

                    # Check for nucleophiles in reactants
                    has_nucleophile = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        or checker.check_fg("Phenol", r)
                        or checker.check_fg("Aromatic thiol", r)
                        or checker.check_fg("Aliphatic thiol", r)
                        for r in reactants
                    )

                    # Check for activating groups in reactants
                    has_activating_group = any(
                        checker.check_fg("Nitro group", r)
                        or checker.check_fg("Nitrile", r)
                        or checker.check_fg("Ester", r)
                        or checker.check_fg("Ketone", r)
                        for r in reactants
                    )

                    # Check if the product has a new C-N, C-O, or C-S bond
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # If we have both an aryl halide and a nucleophile, and potentially an activating group
                        if has_aryl_halide and has_nucleophile:
                            print(f"Found potential nucleophilic aromatic substitution: {rsmi}")
                            print(
                                f"Aryl halide: {has_aryl_halide}, Nucleophile: {has_nucleophile}, Activating group: {has_activating_group}"
                            )
                            has_snar = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_snar
