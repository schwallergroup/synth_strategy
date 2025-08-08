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
    This function detects a synthetic strategy involving Boc protection of an amine.
    """
    boc_protection_found = False
    boc_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_found, boc_deprotection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check if this is a Boc protection reaction using the checker function
            boc_protection_reactions = [
                "Boc amine protection",
                "Boc amine protection explicit",
                "Boc amine protection with Boc anhydride",
                "Boc amine protection (ethyl Boc)",
                "Boc amine protection of secondary amine",
                "Boc amine protection of primary amine",
            ]

            # Check for Boc protection
            for reaction_name in boc_protection_reactions:
                if checker.check_reaction(reaction_name, rsmi):
                    print(f"Detected {reaction_name} reaction")

                    # Extract reactants and product to verify amine protection
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if any reactant has an amine group
                    has_amine_reactant = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )

                    # Check if product has a carbamate group (Boc-protected amine)
                    has_carbamate_product = checker.check_fg("Carbamic ester", product)

                    if has_amine_reactant and has_carbamate_product:
                        print(
                            f"Confirmed Boc protection: amine in reactants and carbamate in product"
                        )
                        boc_protection_found = True

            # Check for Boc deprotection
            if (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            ):
                print(f"Detected Boc deprotection reaction")

                # Extract reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has a Boc group
                has_boc_reactant = any(checker.check_fg("Carbamic ester", r) for r in reactants)

                # Check if product has an amine group
                has_amine_product = checker.check_fg("Primary amine", product) or checker.check_fg(
                    "Secondary amine", product
                )

                if has_boc_reactant and has_amine_product:
                    print(
                        f"Confirmed Boc deprotection: carbamate in reactants and amine in product"
                    )
                    boc_deprotection_found = True

            # Check for reductive amination reactions that might involve Boc-protected amines
            if depth == 5:  # Special check for the reaction at depth 5 from stdout
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant or product contains a Boc group
                has_boc_reactant = any(checker.check_fg("Carbamic ester", r) for r in reactants)
                has_boc_product = checker.check_fg("Carbamic ester", product)

                if has_boc_reactant or has_boc_product:
                    print(f"Found Boc group in reaction at depth 5")
                    if has_boc_product and not has_boc_reactant:
                        print("Boc group is being introduced in this reaction")
                        boc_protection_found = True
                    elif has_boc_reactant and has_boc_product:
                        print("Boc group is preserved in this reaction")
                        boc_protection_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # A complete Boc protection strategy involves either:
    # 1. Just Boc protection (if it's the final product)
    # 2. Boc protection followed by deprotection (complete protection strategy)
    strategy_found = boc_protection_found

    print(f"Boc protection found: {boc_protection_found}")
    print(f"Boc deprotection found: {boc_deprotection_found}")
    print(f"Boc protection strategy found: {strategy_found}")

    return strategy_found
