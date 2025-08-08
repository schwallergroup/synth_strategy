#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    Detects a synthetic strategy involving nosyl protection and deprotection of an amine.
    """
    nosyl_protected = False
    nosyl_deprotected = False

    def dfs_traverse(node):
        nonlocal nosyl_protected, nosyl_deprotected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nosyl deprotection
                # Nosyl deprotection involves removing o-nitrobenzenesulfonyl group from an amine
                for reactant in reactants:
                    if reactant:
                        try:
                            r_mol = Chem.MolFromSmiles(reactant)
                            if r_mol and checker.check_fg("Sulfonamide", reactant):
                                # Check if the sulfonamide has a nitro group attached to the benzene
                                if checker.check_fg("Nitro group", reactant):
                                    # Check if the product is an amine without the nosyl group
                                    p_mol = Chem.MolFromSmiles(product)
                                    if (
                                        p_mol
                                        and (
                                            checker.check_fg("Primary amine", product)
                                            or checker.check_fg("Secondary amine", product)
                                            or checker.check_fg("Tertiary amine", product)
                                        )
                                        and not checker.check_fg("Sulfonamide", product)
                                    ):
                                        nosyl_deprotected = True
                                        print("Detected nosyl deprotection")
                                        print(f"Reactant: {reactant}")
                                        print(f"Product: {product}")
                        except Exception as e:
                            print(f"Error checking deprotection: {e}")
                            continue

                # Alternative check for nosyl deprotection
                if not nosyl_deprotected:
                    has_nosyl_protected = False
                    has_amine_product = False

                    for reactant in reactants:
                        if (
                            reactant
                            and checker.check_fg("Sulfonamide", reactant)
                            and checker.check_fg("Nitro group", reactant)
                        ):
                            has_nosyl_protected = True
                            break

                    if (
                        has_nosyl_protected
                        and (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        )
                        and not checker.check_fg("Sulfonamide", product)
                    ):
                        nosyl_deprotected = True
                        print("Detected nosyl deprotection (alternative method)")
                        print(f"Product: {product}")

                # Check for nosyl protection
                # Check if this is a sulfonamide synthesis reaction
                if checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                ):
                    # Verify that one reactant has a nitro group (nosyl chloride)
                    has_nosyl_cl = False
                    has_amine = False

                    for reactant in reactants:
                        if reactant:
                            try:
                                if checker.check_fg("Nitro group", reactant) and checker.check_fg(
                                    "Sulfonyl halide", reactant
                                ):
                                    has_nosyl_cl = True
                                    print(f"Found nosyl chloride: {reactant}")
                                if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                    "Secondary amine", reactant
                                ):
                                    has_amine = True
                                    print(f"Found amine: {reactant}")
                            except Exception as e:
                                print(f"Error checking protection reactants: {e}")
                                continue

                    # Check if product has sulfonamide with nitro group
                    if has_nosyl_cl and has_amine:
                        try:
                            if checker.check_fg("Sulfonamide", product) and checker.check_fg(
                                "Nitro group", product
                            ):
                                nosyl_protected = True
                                print("Detected nosyl protection")
                                print(f"Product: {product}")
                        except Exception as e:
                            print(f"Error checking protection product: {e}")

                # Alternative check for nosyl protection if reaction checker doesn't work
                if not nosyl_protected:
                    has_nosyl_cl = False
                    has_amine = False
                    nosyl_cl_reactant = None
                    amine_reactant = None

                    for reactant in reactants:
                        if reactant:
                            try:
                                if checker.check_fg("Nitro group", reactant) and checker.check_fg(
                                    "Sulfonyl halide", reactant
                                ):
                                    has_nosyl_cl = True
                                    nosyl_cl_reactant = reactant
                                if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                    "Secondary amine", reactant
                                ):
                                    has_amine = True
                                    amine_reactant = reactant
                            except Exception as e:
                                print(f"Error in alternative protection check: {e}")
                                continue

                    if has_nosyl_cl and has_amine:
                        try:
                            if checker.check_fg("Sulfonamide", product) and checker.check_fg(
                                "Nitro group", product
                            ):
                                print(f"Detected nosyl protection (alternative method)")
                                print(f"Nosyl chloride: {nosyl_cl_reactant}")
                                print(f"Amine: {amine_reactant}")
                                print(f"Product: {product}")
                                nosyl_protected = True
                        except Exception as e:
                            print(f"Error in alternative protection product check: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Protection detected: {nosyl_protected}")
    print(f"Deprotection detected: {nosyl_deprotected}")

    # Return True if both protection and deprotection are detected
    return nosyl_protected and nosyl_deprotected
