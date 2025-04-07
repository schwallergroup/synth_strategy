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
    This function detects a silicon-based protection-deprotection sequence
    for alcohols (specifically TBDPS protection).
    """
    # Track protected molecules and their protection/deprotection status
    protected_molecules = {}

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for protection: R-OH + silyl-Cl → R-O-silyl
            # Look for alcohol in reactants and silyl ether in product
            reactant_has_alcohol = any(
                checker.check_fg("Primary alcohol", r)
                or checker.check_fg("Secondary alcohol", r)
                or checker.check_fg("Tertiary alcohol", r)
                for r in reactants
                if r
            )

            product_has_silyl = checker.check_fg(
                "Silyl protective group", product
            ) or checker.check_fg("TMS ether protective group", product)

            is_protection_reaction = checker.check_reaction(
                "Alcohol protection with silyl ethers", rsmi
            )

            print(
                f"Depth {depth} - Protection check - Alcohol in reactants: {reactant_has_alcohol}, "
                f"Silyl in product: {product_has_silyl}, Is protection reaction: {is_protection_reaction}"
            )

            if is_protection_reaction or (reactant_has_alcohol and product_has_silyl):
                # Find which reactant has the alcohol
                for r in reactants:
                    if (
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                    ):
                        # Store the protected molecule with its depth
                        mol = Chem.MolFromSmiles(product)
                        if mol:
                            protected_molecules[product] = {
                                "protected": True,
                                "depth": depth,
                            }
                            print(f"Silicon protection detected at depth {depth}")
                            break

            # Check for deprotection: R-O-silyl → R-OH
            # Look for silyl ether in reactants and alcohol in product
            reactant_has_silyl = any(
                checker.check_fg("Silyl protective group", r)
                or checker.check_fg("TMS ether protective group", r)
                for r in reactants
                if r
            )

            product_has_alcohol = (
                checker.check_fg("Primary alcohol", product)
                or checker.check_fg("Secondary alcohol", product)
                or checker.check_fg("Tertiary alcohol", product)
            )

            is_deprotection_reaction = (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (double)", rsmi
                )
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (diol)", rsmi
                )
            )

            print(
                f"Depth {depth} - Deprotection check - Silyl in reactants: {reactant_has_silyl}, "
                f"Alcohol in product: {product_has_alcohol}, Is deprotection reaction: {is_deprotection_reaction}"
            )

            if is_deprotection_reaction or (reactant_has_silyl and product_has_alcohol):
                # Find which reactant has the silyl group
                for r in reactants:
                    if checker.check_fg(
                        "Silyl protective group", r
                    ) or checker.check_fg("TMS ether protective group", r):
                        # Store the deprotected molecule with its depth
                        mol = Chem.MolFromSmiles(product)
                        if mol:
                            # Mark as deprotected
                            protected_molecules[r] = protected_molecules.get(r, {})
                            protected_molecules[r]["deprotected"] = True
                            protected_molecules[r]["deprotection_depth"] = depth
                            print(f"Silicon deprotection detected at depth {depth}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least one molecule that was both protected and deprotected
    protection_deprotection_found = False
    for mol, status in protected_molecules.items():
        if status.get("protected") and status.get("deprotected"):
            protection_deprotection_found = True
            print(f"Found molecule that was both protected and deprotected: {mol}")
            break

    # If we didn't find a specific molecule with both steps, check if we found both steps at all
    if not protection_deprotection_found:
        protection_found = any(
            status.get("protected") for status in protected_molecules.values()
        )
        deprotection_found = any(
            status.get("deprotected") for status in protected_molecules.values()
        )
        protection_deprotection_found = protection_found and deprotection_found

    print(
        f"Final result - Protection-deprotection strategy found: {protection_deprotection_found}"
    )

    return protection_deprotection_found
