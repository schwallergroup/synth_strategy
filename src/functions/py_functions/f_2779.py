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
    This function detects a ketone protection/deprotection strategy using ketal.
    It looks for the presence of ketal groups that are later deprotected to ketones.
    """
    ketal_found = False
    ketal_deprotection_found = False
    in_stock_ketal = False

    def dfs_traverse(node, depth=0):
        nonlocal ketal_found, ketal_deprotection_found, in_stock_ketal

        # Check if this is an in-stock molecule with a ketal group
        if node["type"] == "mol" and node.get("in_stock", False):
            if checker.check_fg("Acetal/Ketal", node["smiles"]):
                in_stock_ketal = True
                print(
                    f"In-stock molecule with ketal found at depth {depth}: {node['smiles']}"
                )

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for ketal protection reaction (ketone to ketal)
                    if checker.check_reaction(
                        "Aldehyde or ketone acetalization", rsmi
                    ) or checker.check_reaction("Diol acetalization", rsmi):
                        # Verify that a ketone is being protected
                        reactant_has_ketone = any(
                            checker.check_fg("Ketone", r) for r in reactants if r
                        )
                        product_has_ketal = checker.check_fg("Acetal/Ketal", product)

                        if product_has_ketal:
                            ketal_found = True
                            print(f"Ketal protection found at depth {depth}: {rsmi}")

                    # Check for ketal deprotection reaction (ketal to ketone)
                    if checker.check_reaction(
                        "Ketal hydrolysis to ketone", rsmi
                    ) or checker.check_reaction("Acetal hydrolysis to ketone", rsmi):
                        # Verify that a ketal is being deprotected to a ketone
                        reactant_has_ketal = any(
                            checker.check_fg("Acetal/Ketal", r) for r in reactants if r
                        )
                        product_has_ketone = checker.check_fg("Ketone", product)

                        if reactant_has_ketal and product_has_ketone:
                            ketal_deprotection_found = True
                            print(f"Ketal deprotection found at depth {depth}: {rsmi}")

                    # Fallback method if reaction checker doesn't identify the reactions
                    if not ketal_found:
                        # Check for ketal protection (product has ketal, reactant has ketone or alcohols)
                        reactant_has_ketone = any(
                            checker.check_fg("Ketone", r) for r in reactants if r
                        )
                        reactant_has_alcohols = any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            for r in reactants
                            if r
                        )
                        product_has_ketal = checker.check_fg("Acetal/Ketal", product)

                        if (
                            product_has_ketal
                            and (reactant_has_ketone or reactant_has_alcohols)
                            and not any(
                                checker.check_fg("Acetal/Ketal", r)
                                for r in reactants
                                if r
                            )
                        ):
                            ketal_found = True
                            print(
                                f"Ketal protection found (fallback) at depth {depth}: {rsmi}"
                            )

                    if not ketal_deprotection_found:
                        # Check for ketal deprotection (reactant has ketal, product has ketone)
                        reactant_has_ketal = any(
                            checker.check_fg("Acetal/Ketal", r) for r in reactants if r
                        )
                        product_has_ketone = checker.check_fg("Ketone", product)

                        if (
                            reactant_has_ketal
                            and product_has_ketone
                            and not checker.check_fg("Acetal/Ketal", product)
                        ):
                            ketal_deprotection_found = True
                            print(
                                f"Ketal deprotection found (fallback) at depth {depth}: {rsmi}"
                            )

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if both protection and deprotection are found
    # OR if an in-stock ketal is deprotected (protection happened before synthesis)
    print(f"Ketal protection found: {ketal_found}")
    print(f"In-stock ketal found: {in_stock_ketal}")
    print(f"Ketal deprotection found: {ketal_deprotection_found}")

    return (ketal_found or in_stock_ketal) and ketal_deprotection_found
