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
    This function detects multiple protection-deprotection sequences
    in the synthetic route.
    """
    protection_count = 0
    deprotection_count = 0

    # List of protection reaction types
    protection_reactions = [
        "Alcohol protection with silyl ethers",
        "Protection of carboxylic acid",
        "Boc amine protection",
        "Boc amine protection explicit",
        "Boc amine protection with Boc anhydride",
        "Boc amine protection (ethyl Boc)",
        "Boc amine protection of secondary amine",
        "Boc amine protection of primary amine",
    ]

    # List of deprotection reaction types
    deprotection_reactions = [
        "Alcohol deprotection from silyl ethers",
        "Alcohol deprotection from silyl ethers (double)",
        "Alcohol deprotection from silyl ethers (diol)",
        "Deprotection of carboxylic acid",
        "Boc amine deprotection",
        "Boc amine deprotection of guanidine",
        "Boc amine deprotection to NH-NH2",
        "Ester saponification (methyl deprotection)",
        "Ester saponification (alkyl deprotection)",
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
        "COOH ethyl deprotection",
        "Tert-butyl deprotection of amine",
        "TMS deprotection from alkyne",
        "N-glutarimide deprotection",
        "Phthalimide deprotection",
    ]

    def dfs_traverse(node):
        nonlocal protection_count, deprotection_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for protection reactions
                for protection_type in protection_reactions:
                    if checker.check_reaction(protection_type, rsmi):
                        protection_count += 1
                        print(f"Protection reaction detected: {protection_type} - {rsmi}")
                        break

                # Check for deprotection reactions
                for deprotection_type in deprotection_reactions:
                    if checker.check_reaction(deprotection_type, rsmi):
                        deprotection_count += 1
                        print(f"Deprotection reaction detected: {deprotection_type} - {rsmi}")
                        break

                # Additional check for other protection/deprotection reactions not in the lists
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for silyl protection not caught by reaction checkers
                if any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants
                ) and checker.check_fg("Silyl protective group", product):
                    protection_count += 1
                    print(f"Additional silyl protection detected: {rsmi}")

                # Check for alcohol deprotection not caught by reaction checkers
                if checker.check_fg("Silyl protective group", reactants[0]) and (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                ):
                    deprotection_count += 1
                    print(f"Additional silyl deprotection detected: {rsmi}")

                # Check for Boc protection not caught by reaction checkers
                if (
                    checker.check_fg("Primary amine", reactants[0])
                    or checker.check_fg("Secondary amine", reactants[0])
                ) and checker.check_fg("Boc", product):
                    protection_count += 1
                    print(f"Additional Boc protection detected: {rsmi}")

                # Check for Boc deprotection not caught by reaction checkers
                if checker.check_fg("Boc", reactants[0]) and (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                ):
                    deprotection_count += 1
                    print(f"Additional Boc deprotection detected: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if we have at least 2 protection or deprotection reactions
    has_multiple_protections = (protection_count + deprotection_count) >= 2
    print(f"Protection count: {protection_count}, Deprotection count: {deprotection_count}")
    return has_multiple_protections
