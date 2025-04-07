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
    Detects a strategy involving sequential functionalization of a pyrimidine core,
    particularly through SNAr reactions.
    """
    pyrimidine_snar_count = 0
    has_cyano_group = False

    def dfs_traverse(node, depth=0):
        nonlocal pyrimidine_snar_count, has_cyano_group

        # Check for cyano group in final product (molecule node at depth 0)
        if depth == 0 and node["type"] == "mol":
            mol_smiles = node["smiles"]
            if checker.check_fg("Nitrile", mol_smiles) and checker.check_ring(
                "pyrimidine", mol_smiles
            ):
                print(f"Found cyano group on pyrimidine in final product: {mol_smiles}")
                has_cyano_group = True
            elif checker.check_fg("Nitrile", mol_smiles):
                print(f"Found cyano group in final product (not on pyrimidine): {mol_smiles}")
                has_cyano_group = True

        # Check for SNAr reactions on pyrimidine
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains pyrimidine
            if checker.check_ring("pyrimidine", product):
                # Find pyrimidine-containing reactant
                pyrimidine_reactants = [r for r in reactants if checker.check_ring("pyrimidine", r)]

                if pyrimidine_reactants:
                    # Check for common SNAr reaction types
                    snar_reaction_types = [
                        "heteroaromatic_nuc_sub",
                        "nucl_sub_aromatic_ortho_nitro",
                        "nucl_sub_aromatic_para_nitro",
                        "N-arylation",
                        "Buchwald-Hartwig",
                        "Ullmann-Goldberg Substitution amine",
                        "Ullmann-Goldberg Substitution thiol",
                        "Ullmann-Goldberg Substitution aryl alcohol",
                    ]

                    is_snar = any(checker.check_reaction(rxn, rsmi) for rxn in snar_reaction_types)

                    # If not detected by reaction type, check for characteristic patterns
                    if not is_snar:
                        # Check for halogen on pyrimidine in reactant
                        for pyrimidine_r in pyrimidine_reactants:
                            if (
                                checker.check_fg("Aromatic halide", pyrimidine_r)
                                or checker.check_fg("Primary halide", pyrimidine_r)
                                or checker.check_fg("Secondary halide", pyrimidine_r)
                                or checker.check_fg("Tertiary halide", pyrimidine_r)
                                or checker.check_fg("Alkenyl halide", pyrimidine_r)
                                or checker.check_fg("Triflate", pyrimidine_r)
                                or checker.check_fg("Mesylate", pyrimidine_r)
                                or checker.check_fg("Tosylate", pyrimidine_r)
                            ):

                                # Check if the reaction involves nucleophile addition
                                nucleophiles = [
                                    "Aniline",
                                    "Primary amine",
                                    "Secondary amine",
                                    "Tertiary amine",
                                    "Phenol",
                                    "Primary alcohol",
                                    "Secondary alcohol",
                                    "Tertiary alcohol",
                                    "Aromatic thiol",
                                    "Aliphatic thiol",
                                ]

                                has_nucleophile = any(
                                    any(checker.check_fg(nuc, r) for r in reactants)
                                    for nuc in nucleophiles
                                )

                                if has_nucleophile:
                                    is_snar = True
                                    break

                    if is_snar:
                        print(f"Found pyrimidine SNAr reaction at depth {depth}: {rsmi}")
                        pyrimidine_snar_count += 1

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    print("Starting traversal of synthetic route")
    dfs_traverse(route)

    # Print final counts
    print(f"Total pyrimidine SNAr reactions: {pyrimidine_snar_count}")
    print(f"Has cyano group in final product: {has_cyano_group}")

    # Return True if strategy criteria are met - modified to require only 1 SNAr reaction
    return pyrimidine_snar_count >= 1 and has_cyano_group
