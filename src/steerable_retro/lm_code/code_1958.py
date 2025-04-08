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

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis includes chain extension via nucleophilic substitution.
    """
    found_nucleophilic_substitution = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nucleophilic_substitution

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nucleophilic substitution reactions
            # First, check for specific reaction types
            if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                print(f"Williamson Ether Synthesis detected at depth {depth}")
                found_nucleophilic_substitution = True
            elif checker.check_reaction("S-alkylation of thiols", rsmi) or checker.check_reaction(
                "S-alkylation of thiols (ethyl)", rsmi
            ):
                print(f"S-alkylation of thiols detected at depth {depth}")
                found_nucleophilic_substitution = True
            elif checker.check_reaction(
                "N-alkylation of primary amines with alkyl halides", rsmi
            ) or checker.check_reaction(
                "N-alkylation of secondary amines with alkyl halides", rsmi
            ):
                print(f"N-alkylation with alkyl halides detected at depth {depth}")
                found_nucleophilic_substitution = True
            else:
                # If no specific reaction type is found, check for the presence of reactants and products
                # that would indicate a nucleophilic substitution for chain extension

                # Check for leaving groups in reactants
                has_primary_halide = any(checker.check_fg("Primary halide", r) for r in reactants)
                has_secondary_halide = any(
                    checker.check_fg("Secondary halide", r) for r in reactants
                )
                has_tertiary_halide = any(checker.check_fg("Tertiary halide", r) for r in reactants)
                has_mesylate = any(checker.check_fg("Mesylate", r) for r in reactants)
                has_tosylate = any(checker.check_fg("Tosylate", r) for r in reactants)
                has_triflate = any(checker.check_fg("Triflate", r) for r in reactants)

                # Check for nucleophiles in reactants
                has_alcohol = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Aromatic alcohol", r)
                    for r in reactants
                )
                has_thiol = any(
                    checker.check_fg("Aromatic thiol", r) or checker.check_fg("Aliphatic thiol", r)
                    for r in reactants
                )
                has_amine = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                )

                # Check for products of nucleophilic substitution
                has_ether_product = checker.check_fg("Ether", product)
                has_sulfide_product = checker.check_fg("Monosulfide", product)
                has_amine_product = checker.check_fg("Tertiary amine", product) or checker.check_fg(
                    "Secondary amine", product
                )

                # Verify chain extension via nucleophilic substitution
                if (
                    has_primary_halide
                    or has_secondary_halide
                    or has_tertiary_halide
                    or has_mesylate
                    or has_tosylate
                    or has_triflate
                ):

                    if has_alcohol and has_ether_product:
                        print(
                            f"Chain extension via alcohol nucleophilic substitution detected at depth {depth}"
                        )
                        found_nucleophilic_substitution = True
                    elif has_thiol and has_sulfide_product:
                        print(
                            f"Chain extension via thiol nucleophilic substitution detected at depth {depth}"
                        )
                        found_nucleophilic_substitution = True
                    elif has_amine and has_amine_product:
                        print(
                            f"Chain extension via amine nucleophilic substitution detected at depth {depth}"
                        )
                        found_nucleophilic_substitution = True

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_nucleophilic_substitution
