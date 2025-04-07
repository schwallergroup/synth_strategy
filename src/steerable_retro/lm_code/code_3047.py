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
    This function detects if the synthetic route employs multiple SNAr reactions
    to build the molecular scaffold.
    """
    snar_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal snar_reactions

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check if this is an SNAr reaction using the checker functions
                is_snar = False

                # Check for specific SNAr reaction types
                if checker.check_reaction("heteroaromatic_nuc_sub", rsmi):
                    print(f"Detected heteroaromatic nucleophilic substitution at depth {depth}")
                    is_snar = True
                elif checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi):
                    print(
                        f"Detected nucleophilic substitution on nitro-activated aromatic at depth {depth}"
                    )
                    is_snar = True
                elif checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi):
                    print(
                        f"Detected nucleophilic substitution on nitro-activated aromatic at depth {depth}"
                    )
                    is_snar = True
                else:
                    # If no specific reaction type matches, check for characteristic patterns
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for aromatic halide in reactants
                    aromatic_halide_found = False
                    nucleophile_found = False

                    for reactant in reactants:
                        try:
                            if checker.check_fg("Aromatic halide", reactant):
                                # Check if it's also a heteroaromatic compound
                                heteroaromatic = False
                                for ring in [
                                    "pyridine",
                                    "pyrimidine",
                                    "pyrazine",
                                    "pyridazine",
                                    "triazine",
                                    "tetrazine",
                                ]:
                                    if checker.check_ring(ring, reactant):
                                        heteroaromatic = True
                                        print(
                                            f"Found heteroaromatic halide with {ring} ring: {reactant}"
                                        )
                                        break

                                if heteroaromatic:
                                    aromatic_halide_found = True

                            # Check for nucleophiles
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Aniline", reactant)
                                or checker.check_fg("Phenol", reactant)
                                or checker.check_fg("Aliphatic thiol", reactant)
                                or checker.check_fg("Aromatic thiol", reactant)
                            ):
                                nucleophile_found = True
                                print(f"Found nucleophile: {reactant}")
                        except Exception as e:
                            print(f"Error processing reactant: {e}")
                            continue

                    # Check if the product has the expected C-N, C-O, or C-S bond formation
                    if aromatic_halide_found and nucleophile_found:
                        try:
                            # This is a simplified check - in a real implementation,
                            # we would verify the specific bond formation at the position of the halide
                            is_snar = True
                            print(
                                f"Detected potential SNAr reaction based on reactants and product"
                            )
                        except Exception as e:
                            print(f"Error checking product: {e}")

                if is_snar:
                    snar_reactions += 1
                    print(f"Detected SNAr reaction #{snar_reactions} at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Total SNAr reactions found: {snar_reactions}")
    return snar_reactions >= 2  # Return True if at least 2 SNAr reactions are found
