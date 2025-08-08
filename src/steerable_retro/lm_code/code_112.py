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
    This function detects if the synthesis involves multiple different
    protection/deprotection steps (at least 2 different types).
    """
    protection_types_found = set()

    # Define protection/deprotection reaction types to check
    protection_reactions = [
        # Boc protection/deprotection
        "Boc amine protection",
        "Boc amine protection explicit",
        "Boc amine protection with Boc anhydride",
        "Boc amine protection (ethyl Boc)",
        "Boc amine protection of secondary amine",
        "Boc amine protection of primary amine",
        "Boc amine deprotection",
        "Boc amine deprotection of guanidine",
        "Boc amine deprotection to NH-NH2",
        # Silyl protection/deprotection
        "Alcohol protection with silyl ethers",
        "Alcohol deprotection from silyl ethers",
        "Alcohol deprotection from silyl ethers (double)",
        "Alcohol deprotection from silyl ethers (diol)",
        # Benzyl protection/deprotection
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
        # Other protection/deprotection
        "Cleavage of methoxy ethers to alcohols",
        "Cleavage of alkoxy ethers to alcohols",
        "Ether cleavage to primary alcohol",
        "COOH ethyl deprotection",
        "Tert-butyl deprotection of amine",
        "N-glutarimide deprotection",
        "Phthalimide deprotection",
        "TMS deprotection from alkyne",
        "Protection of carboxylic acid",
        "Deprotection of carboxylic acid",
        # Additional protection/deprotection reactions
        "Ester saponification (methyl deprotection)",
        "Ester saponification (alkyl deprotection)",
    ]

    # Group reactions by protection type
    protection_groups = {
        "boc": [
            "Boc amine protection",
            "Boc amine protection explicit",
            "Boc amine protection with Boc anhydride",
            "Boc amine protection (ethyl Boc)",
            "Boc amine protection of secondary amine",
            "Boc amine protection of primary amine",
            "Boc amine deprotection",
            "Boc amine deprotection of guanidine",
            "Boc amine deprotection to NH-NH2",
        ],
        "silyl": [
            "Alcohol protection with silyl ethers",
            "Alcohol deprotection from silyl ethers",
            "Alcohol deprotection from silyl ethers (double)",
            "Alcohol deprotection from silyl ethers (diol)",
        ],
        "benzyl": ["Hydroxyl benzyl deprotection", "Carboxyl benzyl deprotection"],
        "ether": [
            "Cleavage of methoxy ethers to alcohols",
            "Cleavage of alkoxy ethers to alcohols",
            "Ether cleavage to primary alcohol",
        ],
        "carboxyl": [
            "COOH ethyl deprotection",
            "Protection of carboxylic acid",
            "Deprotection of carboxylic acid",
            "Ester saponification (methyl deprotection)",
            "Ester saponification (alkyl deprotection)",
        ],
        "amine": [
            "Tert-butyl deprotection of amine",
            "N-glutarimide deprotection",
            "Phthalimide deprotection",
        ],
        "alkyne": ["TMS deprotection from alkyne"],
    }

    # Reverse mapping from reaction to group
    reaction_to_group = {}
    for group, reactions in protection_groups.items():
        for reaction in reactions:
            reaction_to_group[reaction] = group

    # Additional check for protection/deprotection by functional group changes
    def check_protection_by_fg(rsmi):
        try:
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            # Check for Boc group
            if (
                checker.check_fg("Boc", reactants_part)
                and not checker.check_fg("Boc", products_part)
            ) or (
                not checker.check_fg("Boc", reactants_part)
                and checker.check_fg("Boc", products_part)
            ):
                return "boc"

            # Check for silyl groups
            if (
                (
                    checker.check_fg("Silyl protective group", reactants_part)
                    and not checker.check_fg("Silyl protective group", products_part)
                )
                or (
                    not checker.check_fg("Silyl protective group", reactants_part)
                    and checker.check_fg("Silyl protective group", products_part)
                )
                or (
                    checker.check_fg("TMS ether protective group", reactants_part)
                    and not checker.check_fg("TMS ether protective group", products_part)
                )
                or (
                    not checker.check_fg("TMS ether protective group", reactants_part)
                    and checker.check_fg("TMS ether protective group", products_part)
                )
            ):
                return "silyl"

            # Check for benzyl protection/deprotection
            if (
                checker.check_fg("Ether", reactants_part)
                and checker.check_fg("Phenol", products_part)
                and not checker.check_fg("Phenol", reactants_part)
            ):
                # Check if benzyl group is involved
                if "OCH2Ph" in reactants_part or "OCH2c1ccccc1" in reactants_part:
                    return "benzyl"

            # Check for ester deprotection (carboxyl protection)
            if (
                checker.check_fg("Ester", reactants_part)
                and checker.check_fg("Carboxylic acid", products_part)
                and not checker.check_fg("Carboxylic acid", reactants_part)
            ):
                return "carboxyl"

            # Check for methyl/ethyl ether deprotection
            if (
                checker.check_fg("Ether", reactants_part)
                and checker.check_fg("Primary alcohol", products_part)
                and not checker.check_fg("Primary alcohol", reactants_part)
            ):
                return "ether"

            # Check for amine protection/deprotection
            if (
                checker.check_fg("Primary amine", reactants_part)
                and not checker.check_fg("Primary amine", products_part)
            ) or (
                not checker.check_fg("Primary amine", reactants_part)
                and checker.check_fg("Primary amine", products_part)
            ):
                # Check if not already counted as Boc
                if not checker.check_fg("Boc", reactants_part) and not checker.check_fg(
                    "Boc", products_part
                ):
                    return "amine"

            # Check for TMS alkyne protection/deprotection
            if (
                checker.check_fg("Alkyne", reactants_part)
                and not checker.check_fg("Alkyne", products_part)
            ) or (
                not checker.check_fg("Alkyne", reactants_part)
                and checker.check_fg("Alkyne", products_part)
            ):
                if (
                    "TMS" in reactants_part
                    or "Si(CH3)3" in reactants_part
                    or "TMS" in products_part
                    or "Si(CH3)3" in products_part
                ):
                    return "alkyne"

            return None
        except Exception as e:
            print(f"Error in check_protection_by_fg: {e}")
            return None

    def check_reaction_by_group(rsmi):
        """Check for protection/deprotection reactions by group"""
        for group, reactions in protection_groups.items():
            for reaction in reactions:
                if checker.check_reaction(reaction, rsmi):
                    return group
        return None

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi[:50]}...")

            # First check by reaction type (grouped for efficiency)
            group = check_reaction_by_group(rsmi)
            if group:
                protection_types_found.add(group)
                print(f"Found protection/deprotection reaction group: {group}")
            else:
                # If no match found by reaction type, try checking by functional group changes
                group = check_protection_by_fg(rsmi)
                if group:
                    protection_types_found.add(group)
                    print(f"Found protection/deprotection by functional group analysis: {group}")

            # Additional check for specific reaction types not covered by groups
            # Check for ester hydrolysis (saponification)
            if (
                not group
                and checker.check_fg("Ester", rsmi.split(">")[0])
                and checker.check_fg("Carboxylic acid", rsmi.split(">")[-1])
            ):
                protection_types_found.add("carboxyl")
                print(f"Found ester hydrolysis (carboxyl deprotection)")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check for specific patterns in the entire route
    def check_entire_route(node, patterns_found=None):
        if patterns_found is None:
            patterns_found = set()

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for protected molecules
            if checker.check_fg("Boc", mol_smiles):
                patterns_found.add("boc")
            if checker.check_fg("Silyl protective group", mol_smiles) or checker.check_fg(
                "TMS ether protective group", mol_smiles
            ):
                patterns_found.add("silyl")
            if checker.check_fg("Ester", mol_smiles) and not node.get("in_stock", False):
                # Only count non-starting material esters as potential protection
                patterns_found.add("carboxyl")

        # Recursively check children
        for child in node.get("children", []):
            check_entire_route(child, patterns_found)

        return patterns_found

    # Add any protection patterns found in the entire route
    additional_patterns = check_entire_route(route)
    protection_types_found.update(additional_patterns)

    result = len(protection_types_found) >= 2
    print(
        f"Multiple protection/deprotection types detected: {result} - Types: {protection_types_found}"
    )
    return result
