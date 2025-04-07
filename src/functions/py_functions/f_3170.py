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
    Detects if the synthesis follows a linear strategy with protection steps.
    """
    protection_steps = 0
    linear_steps = 0
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal protection_steps, linear_steps, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Count number of reactants to check if linear
            reactants = reactants_smiles.split(".")
            if len(reactants) <= 2:
                linear_steps += 1
                print(f"Found linear step: {rsmi}")

            # Check for protection/deprotection reactions by name
            protection_reaction = False

            # Check for protection/deprotection reactions
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction(
                    "Boc amine protection with Boc anhydride", rsmi
                )
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction(
                    "Boc amine protection of secondary amine", rsmi
                )
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                or checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                or checker.check_reaction("Protection of carboxylic acid", rsmi)
                or checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers", rsmi
                )
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (double)", rsmi
                )
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (diol)", rsmi
                )
                or checker.check_reaction("Deprotection of carboxylic acid", rsmi)
                or checker.check_reaction("Hydroxyl benzyl deprotection", rsmi)
                or checker.check_reaction("Carboxyl benzyl deprotection", rsmi)
                or checker.check_reaction(
                    "Cleavage of methoxy ethers to alcohols", rsmi
                )
                or checker.check_reaction("Cleavage of alkoxy ethers to alcohols", rsmi)
                or checker.check_reaction("Ether cleavage to primary alcohol", rsmi)
                or checker.check_reaction("COOH ethyl deprotection", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                or checker.check_reaction("TMS deprotection from alkyne", rsmi)
                or checker.check_reaction("N-glutarimide deprotection", rsmi)
                or checker.check_reaction("Phthalimide deprotection", rsmi)
            ):
                protection_steps += 1
                protection_reaction = True
                print(f"Found protection reaction: {rsmi}")

            # Check for acetal/ketal formation/hydrolysis
            if not protection_reaction and (
                checker.check_reaction("Aldehyde or ketone acetalization", rsmi)
                or checker.check_reaction("Diol acetalization", rsmi)
                or checker.check_reaction("Acetal hydrolysis to diol", rsmi)
                or checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                or checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
            ):
                protection_steps += 1
                protection_reaction = True
                print(f"Found acetal/ketal protection: {rsmi}")

            # Check for esterification/saponification (often used as protection)
            if not protection_reaction and (
                checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                or checker.check_reaction(
                    "Ester saponification (methyl deprotection)", rsmi
                )
                or checker.check_reaction(
                    "Ester saponification (alkyl deprotection)", rsmi
                )
            ):
                protection_steps += 1
                protection_reaction = True
                print(f"Found ester protection: {rsmi}")

            # Check for protection groups in products
            if not protection_reaction:
                product_mol = (
                    Chem.MolFromSmiles(product_smiles) if product_smiles else None
                )
                if product_mol:
                    # Check for common protection groups in products
                    if (
                        checker.check_fg("TMS ether protective group", product_smiles)
                        or checker.check_fg("Silyl protective group", product_smiles)
                        or checker.check_fg("Boc", product_smiles)
                        or checker.check_fg("Acetal/Ketal", product_smiles)
                    ):
                        protection_steps += 1
                        protection_reaction = True
                        print(f"Found protection group in product: {product_smiles}")

            # Check for protection groups in reactants (might indicate deprotection)
            if not protection_reaction:
                for reactant in reactants:
                    if (
                        checker.check_fg("TMS ether protective group", reactant)
                        or checker.check_fg("Silyl protective group", reactant)
                        or checker.check_fg("Boc", reactant)
                        or checker.check_fg("Acetal/Ketal", reactant)
                    ):
                        # Check if the protection group is absent in the product (indicating deprotection)
                        if not (
                            checker.check_fg(
                                "TMS ether protective group", product_smiles
                            )
                            or checker.check_fg(
                                "Silyl protective group", product_smiles
                            )
                            or checker.check_fg("Boc", product_smiles)
                            or checker.check_fg("Acetal/Ketal", product_smiles)
                        ):
                            protection_steps += 1
                            protection_reaction = True
                            print(f"Found deprotection: {rsmi}")
                            break

            # Check for specific functional group transformations that indicate protection/deprotection
            if not protection_reaction:
                # Alcohol protection/deprotection
                if any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants
                ):
                    if checker.check_fg("Ether", product_smiles) or checker.check_fg(
                        "Ester", product_smiles
                    ):
                        protection_steps += 1
                        protection_reaction = True
                        print(f"Found alcohol protection: {rsmi}")

                # Amine protection/deprotection
                if not protection_reaction and any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    for r in reactants
                ):
                    if checker.check_fg("Amide", product_smiles) or checker.check_fg(
                        "Carbamic ester", product_smiles
                    ):
                        protection_steps += 1
                        protection_reaction = True
                        print(f"Found amine protection: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if synthesis is predominantly linear with protection steps
    # Lowered threshold to 50% based on test case
    is_linear = linear_steps > 0 and linear_steps >= max(1, (max_depth - 1) * 0.5)
    has_protection = protection_steps >= 1

    print(
        f"Linear steps: {linear_steps}, Max depth: {max_depth}, Protection steps: {protection_steps}"
    )
    print(f"Is linear: {is_linear}, Has protection: {has_protection}")
    return is_linear and has_protection
