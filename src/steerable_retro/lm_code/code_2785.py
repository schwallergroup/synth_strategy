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
    Detects a strategy involving the conversion of an aldehyde to an alkene.
    Common reactions include Wittig, Julia olefination, and aldol condensation.
    """
    conversion_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal conversion_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Check for aldehyde in reactants
            aldehyde_in_reactants = any(checker.check_fg("Aldehyde", r) for r in reactants if r)

            # Check for alkene-related functional groups in product
            alkene_in_product = False
            if product:
                alkene_in_product = (
                    checker.check_fg("Alkene", product)
                    or checker.check_fg("Vinyl", product)
                    or checker.check_fg("Allyl", product)
                    or checker.check_fg("Conjugated diene", product)
                )

            # Check for specific reactions that convert aldehydes to alkenes
            wittig_reaction = (
                checker.check_reaction("Wittig reaction", rsmi)
                or checker.check_reaction("Wittig", rsmi)
                or checker.check_reaction("{Wittig}", rsmi)
            )
            julia_reaction = checker.check_reaction("Julia Olefination", rsmi)
            aldol_condensation = checker.check_reaction("Aldol condensation", rsmi)

            print(f"  Aldehyde in reactants: {aldehyde_in_reactants}")
            print(f"  Alkene in product: {alkene_in_product}")
            print(f"  Wittig reaction: {wittig_reaction}")
            print(f"  Julia reaction: {julia_reaction}")
            print(f"  Aldol condensation: {aldol_condensation}")

            # If we have an aldehyde in reactants and an alkene in product, this is likely our target conversion
            if aldehyde_in_reactants and alkene_in_product:
                print(f"  Found aldehyde in reactants and alkene in product")

                # Check for specific named reactions
                if wittig_reaction:
                    print("  Detected Wittig reaction converting aldehyde to alkene")
                    conversion_detected = True
                elif julia_reaction:
                    print("  Detected Julia olefination converting aldehyde to alkene")
                    conversion_detected = True
                elif aldol_condensation:
                    print("  Detected aldol condensation converting aldehyde to alkene")
                    conversion_detected = True
                elif checker.check_reaction(
                    "Aldehyde or ketone to alpha,beta-unsaturated carbonyl", rsmi
                ):
                    print("  Detected conversion of aldehyde to alpha,beta-unsaturated carbonyl")
                    conversion_detected = True
                else:
                    # Check for other reactions or patterns
                    # Look for patterns in atom-mapped SMILES that indicate aldehyde to alkene conversion
                    try:
                        # Find aldehyde reactant
                        aldehyde_reactant = None
                        for r in reactants:
                            if checker.check_fg("Aldehyde", r):
                                aldehyde_reactant = r
                                break

                        if aldehyde_reactant:
                            # Check if this is a reaction that typically converts aldehydes to alkenes
                            # Look for patterns in atom-mapped SMILES
                            if "=" in product and "[CH" in product:
                                print(
                                    "  Detected potential aldehyde to alkene conversion through unspecified reaction"
                                )
                                conversion_detected = True
                    except Exception as e:
                        print(f"  Error in pattern check: {e}")

            # Special case for Wittig-like reactions where the product might not be detected as an alkene
            if aldehyde_in_reactants and not conversion_detected:
                if "CH2=" in product or "=[CH]" in product or "=[CH2]" in product:
                    print("  Detected potential aldehyde to alkene conversion (special pattern)")
                    conversion_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: Aldehyde to alkene conversion detected: {conversion_detected}")
    return conversion_detected
