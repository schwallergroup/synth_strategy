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
    This function detects a strategy involving olefination followed by C-N bond disconnection.
    In retrosynthetic analysis, this means C-N disconnection appears at a lower depth than olefination.
    """
    has_olefination = False
    has_cn_disconnection = False
    olefination_depth = -1
    cn_disconnection_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_olefination, has_cn_disconnection, olefination_depth, cn_disconnection_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for olefination reactions
                if not has_olefination:
                    if (
                        checker.check_reaction("Wittig", rsmi)
                        or checker.check_reaction("Wittig with Phosphonium", rsmi)
                        or checker.check_reaction("Julia Olefination", rsmi)
                        or checker.check_reaction("Heck terminal vinyl", rsmi)
                        or checker.check_reaction("Heck_terminal_vinyl", rsmi)
                        or checker.check_reaction("Heck_non-terminal_vinyl", rsmi)
                        or checker.check_reaction(
                            "Olefination of ketones with Grignard reagents", rsmi
                        )
                        or checker.check_reaction(
                            "Olefination of aldehydes with Grignard reagents", rsmi
                        )
                    ):
                        has_olefination = True
                        olefination_depth = depth
                        print(f"Detected olefination reaction at depth {depth}: {rsmi}")

                    # Fallback to checking for C=C bond formation if specific reaction not found
                    elif not has_olefination:
                        try:
                            product_mol = Chem.MolFromSmiles(product)

                            # Check all reactants
                            for reactant in reactants:
                                reactant_mol = Chem.MolFromSmiles(reactant)

                                if reactant_mol is not None and product_mol is not None:
                                    # Count C=C bonds in reactant and product
                                    cc_double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
                                    reactant_cc_bonds = len(
                                        reactant_mol.GetSubstructMatches(cc_double_bond_pattern)
                                    )
                                    product_cc_bonds = len(
                                        product_mol.GetSubstructMatches(cc_double_bond_pattern)
                                    )

                                    if product_cc_bonds > reactant_cc_bonds:
                                        has_olefination = True
                                        olefination_depth = depth
                                        print(
                                            f"Detected olefination by C=C bond increase at depth {depth}"
                                        )
                                        break
                        except Exception as e:
                            print(f"Error in olefination detection: {e}")

                # Check for C-N bond disconnection (independent of olefination)
                if not has_cn_disconnection:
                    # Check for common C-N disconnection reactions
                    cn_disconnection_reactions = [
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        "Hydrolysis of amides/imides/carbamates",
                        "Hydrogenolysis of amides/imides/carbamates",
                        "Boc amine deprotection",
                        "N-glutarimide deprotection",
                        "Phthalimide deprotection",
                        "Reductive amination with aldehyde",
                        "Reductive amination with ketone",
                        "Reductive amination with alcohol",
                        "Reductive amination",
                        "N-alkylation of primary amines with alkyl halides",
                        "N-alkylation of secondary amines with alkyl halides",
                        "Alkylation of amines",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                        "Buchwald-Hartwig",
                        "Urea synthesis via isocyanate and primary amine",
                        "Urea synthesis via isocyanate and secondary amine",
                    ]

                    for rxn_type in cn_disconnection_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            has_cn_disconnection = True
                            cn_disconnection_depth = depth
                            print(
                                f"Detected C-N bond disconnection reaction at depth {depth}: {rxn_type}"
                            )
                            break

                    # If no specific reaction found, check for C-N bond breaking
                    if not has_cn_disconnection:
                        try:
                            # Check if product has nitrogen
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and "N" in product:
                                # Check for C-N bond formation in forward direction (disconnection in retro)
                                cn_bond_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                                product_cn_bonds = len(
                                    product_mol.GetSubstructMatches(cn_bond_pattern)
                                )

                                # Check if any reactant has fewer C-N bonds
                                for reactant in reactants:
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol is None:
                                        continue

                                    reactant_cn_bonds = len(
                                        reactant_mol.GetSubstructMatches(cn_bond_pattern)
                                    )

                                    # If product has more C-N bonds than reactant, it's a C-N bond formation
                                    if product_cn_bonds > reactant_cn_bonds:
                                        has_cn_disconnection = True
                                        cn_disconnection_depth = depth
                                        print(
                                            f"Detected C-N bond formation (disconnection in retro) at depth {depth}"
                                        )
                                        break

                                    # Also check for nitrogen-containing functional groups
                                    n_containing_fgs = [
                                        "Primary amine",
                                        "Secondary amine",
                                        "Tertiary amine",
                                        "Primary amide",
                                        "Secondary amide",
                                        "Tertiary amide",
                                        "Urea",
                                        "Thiourea",
                                        "Aniline",
                                        "Azide",
                                        "Nitrile",
                                        "Nitro group",
                                        "Isocyanate",
                                        "Isocyanide",
                                        "Sulfonamide",
                                    ]

                                    for fg in n_containing_fgs:
                                        if checker.check_fg(fg, product) and not checker.check_fg(
                                            fg, reactant
                                        ):
                                            has_cn_disconnection = True
                                            cn_disconnection_depth = depth
                                            print(f"Detected {fg} formation at depth {depth}")
                                            break

                                    if has_cn_disconnection:
                                        break
                        except Exception as e:
                            print(f"Error in C-N disconnection detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both reactions were found and in the correct sequence
    # In retrosynthesis, C-N disconnection should be at a lower or equal depth (later or same stage)
    # than olefination (earlier stage)
    correct_sequence = (
        has_olefination and has_cn_disconnection and cn_disconnection_depth <= olefination_depth
    )

    print(f"Olefination found: {has_olefination} at depth {olefination_depth}")
    print(f"C-N disconnection found: {has_cn_disconnection} at depth {cn_disconnection_depth}")
    print(f"Correct sequence (C-N after or with olefination in synthesis): {correct_sequence}")

    return correct_sequence
