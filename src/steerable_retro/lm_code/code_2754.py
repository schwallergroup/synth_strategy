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
    Detects if the synthetic route involves protecting group strategy,
    specifically TMS protection/deprotection for alkynes and Boc protection/deprotection for amines.
    """
    tms_protection = False
    tms_deprotection = False
    boc_protection = False
    boc_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal tms_protection, tms_deprotection, boc_protection, boc_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            reactants = reactants_str.split(".")
            product_str = rsmi.split(">")[-1]

            # In retrosynthesis, the product of the reaction is what we're trying to make
            # and the reactants are the precursors

            # Check for TMS protection/deprotection
            # TMS deprotection: TMS-alkyne → terminal alkyne
            if checker.check_reaction("TMS deprotection from alkyne", rsmi):
                print(f"Detected TMS deprotection reaction at depth {depth}: {rsmi}")
                tms_deprotection = True

            # Check if any reactant has TMS group and product has terminal alkyne
            for reactant in reactants:
                if checker.check_fg("TMS ether protective group", reactant) and checker.check_fg(
                    "Alkyne", product_str
                ):
                    if not checker.check_fg("TMS ether protective group", product_str):
                        print(
                            f"Detected TMS deprotection at depth {depth}: {reactant} -> {product_str}"
                        )
                        tms_deprotection = True

            # Check for silyl protection reactions that might apply to alkynes
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                if checker.check_fg("Alkyne", reactants_str) and checker.check_fg(
                    "Silyl protective group", product_str
                ):
                    print(
                        f"Detected silyl protection at depth {depth}: {reactants_str} -> {product_str}"
                    )
                    tms_protection = True
                elif checker.check_fg("Alkyne", reactants_str) and checker.check_fg(
                    "TMS ether protective group", product_str
                ):
                    print(
                        f"Detected TMS protection at depth {depth}: {reactants_str} -> {product_str}"
                    )
                    tms_protection = True

            # Check if any reactant has terminal alkyne and product has TMS group
            for reactant in reactants:
                if checker.check_fg("Alkyne", reactant) and not checker.check_fg(
                    "TMS ether protective group", reactant
                ):
                    if checker.check_fg(
                        "TMS ether protective group", product_str
                    ) and checker.check_fg("Alkyne", product_str):
                        print(
                            f"Detected TMS protection at depth {depth}: {reactant} -> {product_str}"
                        )
                        tms_protection = True
                    elif checker.check_fg(
                        "Silyl protective group", product_str
                    ) and checker.check_fg("Alkyne", product_str):
                        print(
                            f"Detected silyl protection at depth {depth}: {reactant} -> {product_str}"
                        )
                        tms_protection = True

            # Check for alcohol deprotection from silyl ethers
            if (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
            ):
                print(f"Detected silyl deprotection reaction at depth {depth}: {rsmi}")
                tms_deprotection = True

            # Boc deprotection: Boc-amine → amine
            if (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            ):
                print(f"Detected Boc deprotection reaction at depth {depth}: {rsmi}")
                boc_deprotection = True

            # Check if any reactant has Boc group and product has free amine
            for reactant in reactants:
                if checker.check_fg("Boc", reactant):
                    if checker.check_fg("Primary amine", product_str) or checker.check_fg(
                        "Secondary amine", product_str
                    ):
                        if not checker.check_fg("Boc", product_str):
                            print(
                                f"Detected Boc deprotection at depth {depth}: {reactant} -> {product_str}"
                            )
                            boc_deprotection = True

            # Boc protection: amine → Boc-amine
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                print(f"Detected Boc protection reaction at depth {depth}: {rsmi}")
                boc_protection = True

            # Check if any reactant has free amine and product has Boc group
            for reactant in reactants:
                if (
                    checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                ) and not checker.check_fg("Boc", reactant):
                    if checker.check_fg("Boc", product_str):
                        print(
                            f"Detected Boc protection at depth {depth}: {reactant} -> {product_str}"
                        )
                        boc_protection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return true if we have either protection or deprotection for either TMS or Boc
    tms_strategy = tms_protection or tms_deprotection
    boc_strategy = boc_protection or boc_deprotection

    print(f"TMS strategy detected: {tms_strategy}")
    print(f"Boc strategy detected: {boc_strategy}")

    # Return true if either TMS or Boc protection strategy is detected
    return tms_strategy or boc_strategy
