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
    This function detects a synthetic strategy involving the installation of an alkoxy chain
    linker between functional groups.
    """
    has_alkoxy_linker = False
    potential_linker_molecules = (
        set()
    )  # Track molecules that might contain linker precursors

    def is_alkoxy_chain(mol_smiles):
        """Check if the molecule contains an alkoxy chain of sufficient length"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Look for patterns like -O-C-C-O- or -O-C-C-C-O- (PEG-like)
        peg_like_pattern = Chem.MolFromSmarts("*-O-[CH2]-[CH2]-[CH2]*")
        if mol.HasSubstructMatch(peg_like_pattern):
            return True

        # Look for patterns like -O-C-C-C- (alkoxy chain)
        alkoxy_pattern = Chem.MolFromSmarts("*-O-[CH2]-[CH2]-[CH2]*")
        if mol.HasSubstructMatch(alkoxy_pattern):
            return True

        # Look for shorter alkoxy chains too
        short_alkoxy = Chem.MolFromSmarts("*-O-[CH2]-[CH2]*")
        if mol.HasSubstructMatch(short_alkoxy):
            return True

        return False

    def connects_functional_groups(mol_smiles):
        """Check if the molecule has an alkoxy chain connecting different functional groups"""
        # List of functional groups that might be connected by linkers
        fg_list = [
            "Carboxylic acid",
            "Ester",
            "Amide",
            "Primary amide",
            "Secondary amide",
            "Tertiary amide",
            "Primary amine",
            "Secondary amine",
            "Tertiary amine",
            "Nitro group",
            "Nitrile",
            "Aromatic halide",
            "Phenol",
            "Primary alcohol",
            "Secondary alcohol",
            "Tertiary alcohol",
            "Sulfonamide",
            "Sulfone",
            "Sulfoxide",
        ]

        # Count functional groups
        fg_count = sum(1 for fg in fg_list if checker.check_fg(fg, mol_smiles))

        # If we have at least two functional groups and an alkoxy chain, it's likely a linker
        return fg_count >= 2 and is_alkoxy_chain(mol_smiles)

    def dfs_traverse(node, depth=0):
        nonlocal has_alkoxy_linker

        if node["type"] == "mol":
            # Track molecules that might contain linkers
            if checker.check_fg("Ether", node["smiles"]) and connects_functional_groups(
                node["smiles"]
            ):
                potential_linker_molecules.add(node["smiles"])
                print(f"Found potential molecule with alkoxy linker: {node['smiles']}")

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction creates a new ether bond
            product_has_ether = checker.check_fg("Ether", product)
            all_reactants_have_ether = all(
                checker.check_fg("Ether", r) for r in reactants if r
            )
            new_ether_formed = product_has_ether and not all_reactants_have_ether

            # Check for Williamson ether synthesis and related reactions
            if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                # Check if reactants contain alcohol/phenol and alkyl halide
                has_alcohol = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Phenol", r)
                    for r in reactants
                )

                has_alkyl_halide = any(
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                    for r in reactants
                )

                # Check if the product has a linker-like structure
                if (
                    has_alcohol
                    and has_alkyl_halide
                    and connects_functional_groups(product)
                ):
                    print(
                        f"Found alkoxy chain linker installation via Williamson ether synthesis: {rsmi}"
                    )
                    has_alkoxy_linker = True

            # Check for other ether formation reactions that could install linkers
            elif any(
                checker.check_reaction(rxn_type, rsmi)
                for rxn_type in [
                    "Mitsunobu aryl ether",
                    "Ullmann-Goldberg Substitution aryl alcohol",
                    "Chan-Lam etherification",
                    "Alcohol to ether",
                    "{Williamson ether}",
                ]
            ):
                # Check if product has an ether that wasn't in all reactants
                if new_ether_formed and connects_functional_groups(product):
                    print(
                        f"Found alkoxy chain linker installation via alternative ether formation: {rsmi}"
                    )
                    has_alkoxy_linker = True

            # Check for alkylation of alcohols with diazo compounds
            elif checker.check_reaction(
                "O-alkylation of carboxylic acids with diazo compounds", rsmi
            ) or checker.check_reaction(
                "O-alkylation of amides with diazo compounds", rsmi
            ):
                if new_ether_formed and connects_functional_groups(product):
                    print(
                        f"Found alkoxy chain linker installation via O-alkylation with diazo compounds: {rsmi}"
                    )
                    has_alkoxy_linker = True

            # Check for nucleophilic substitution reactions that could install alkoxy linkers
            elif any(
                checker.check_reaction(rxn_type, rsmi)
                for rxn_type in [
                    "S-alkylation of thiols with alcohols",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Ring opening of epoxide with amine",
                ]
            ):
                if new_ether_formed and connects_functional_groups(product):
                    print(
                        f"Found alkoxy chain linker installation via nucleophilic substitution: {rsmi}"
                    )
                    has_alkoxy_linker = True

            # Check for esterification reactions that could be part of linker installation
            elif any(
                checker.check_reaction(rxn_type, rsmi)
                for rxn_type in [
                    "Esterification of Carboxylic Acids",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Transesterification",
                ]
            ):
                # Check if the product contains both ester and ether groups
                if checker.check_fg("Ester", product) and checker.check_fg(
                    "Ether", product
                ):
                    if connects_functional_groups(product):
                        print(
                            f"Found potential alkoxy chain linker installation via esterification: {rsmi}"
                        )
                        has_alkoxy_linker = True

            # Check for any reaction that produces a product with a linker-like structure
            elif new_ether_formed:
                # Check if the product has a structure consistent with an alkoxy linker
                if connects_functional_groups(product):
                    print(
                        f"Found potential alkoxy chain linker installation via other reaction: {rsmi}"
                    )
                    has_alkoxy_linker = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_alkoxy_linker
