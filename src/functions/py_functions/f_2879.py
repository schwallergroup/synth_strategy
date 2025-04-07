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
    This function detects if the synthetic route involves TBDMS protection of an alcohol.
    """
    tbdms_protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal tbdms_protection_found

        if tbdms_protection_found:
            return  # Early return if already found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for protection or deprotection reactions
                is_protection_reaction = checker.check_reaction(
                    "Alcohol protection with silyl ethers", rsmi
                )
                is_deprotection_reaction = checker.check_reaction(
                    "Alcohol deprotection from silyl ethers", rsmi
                )

                # In retrosynthesis, we might encounter either protection or deprotection
                if is_protection_reaction or is_deprotection_reaction:
                    # Check if it's specifically TBDMS
                    # For protection: TBDMS should be in product
                    # For deprotection: TBDMS should be in reactants
                    mol_to_check = (
                        product if is_protection_reaction else ".".join(reactants)
                    )

                    # Check for silyl protective group
                    has_silyl_group = checker.check_fg(
                        "Silyl protective group", mol_to_check
                    )

                    # Check specifically for TBDMS pattern
                    has_tbdms = False
                    if has_silyl_group:
                        mol = Chem.MolFromSmiles(mol_to_check)
                        if mol:
                            # TBDMS pattern: Si(C)(C)C(C)(C)C
                            tbdms_pattern = Chem.MolFromSmarts("[Si](C)(C)C(C)(C)C")
                            has_tbdms = mol.HasSubstructMatch(tbdms_pattern)
                            print(f"TBDMS pattern found: {has_tbdms}")

                    if has_tbdms:
                        print(
                            f"TBDMS {'protection' if is_protection_reaction else 'deprotection'} detected in reaction: {rsmi}"
                        )
                        tbdms_protection_found = True
                        return

                # Alternative detection method
                # Check for alcohol and silyl group
                has_alcohol_in_reactants = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Aromatic alcohol", r)
                    for r in reactants
                )

                has_alcohol_in_product = (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                    or checker.check_fg("Aromatic alcohol", product)
                )

                has_silyl_in_reactants = any(
                    checker.check_fg("Silyl protective group", r) for r in reactants
                )
                has_silyl_in_product = checker.check_fg(
                    "Silyl protective group", product
                )

                # Check for TBDMS pattern
                has_tbdms_in_reactants = False
                has_tbdms_in_product = False

                if has_silyl_in_reactants:
                    for r in reactants:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            tbdms_pattern = Chem.MolFromSmarts("[Si](C)(C)C(C)(C)C")
                            if mol.HasSubstructMatch(tbdms_pattern):
                                has_tbdms_in_reactants = True
                                break

                if has_silyl_in_product:
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        tbdms_pattern = Chem.MolFromSmarts("[Si](C)(C)C(C)(C)C")
                        has_tbdms_in_product = mol.HasSubstructMatch(tbdms_pattern)

                # In forward direction: alcohol + reagent -> silyl ether (protection)
                # In retrosynthesis: silyl ether -> alcohol + reagent (deprotection)
                if (has_alcohol_in_reactants and has_tbdms_in_product) or (
                    has_tbdms_in_reactants and has_alcohol_in_product
                ):
                    print(
                        f"TBDMS protection/deprotection detected through functional group analysis: {rsmi}"
                    )
                    tbdms_protection_found = True
                    return

                if not tbdms_protection_found:
                    if not is_protection_reaction and not is_deprotection_reaction:
                        print("Not a silyl ether protection/deprotection reaction")
                    if has_silyl_in_product and not has_tbdms_in_product:
                        print(
                            "Silyl group detected in product but not TBDMS specifically"
                        )
                    if has_silyl_in_reactants and not has_tbdms_in_reactants:
                        print(
                            "Silyl group detected in reactants but not TBDMS specifically"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"TBDMS protection strategy found: {tbdms_protection_found}")
    return tbdms_protection_found
