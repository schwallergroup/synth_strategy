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
    Detects if the synthetic route employs a protection-deprotection strategy
    with at least two different protecting groups.
    """
    # Track protection and deprotection reactions
    protection_reactions = {
        "boc": False,
        "silyl": False,
        "benzyl": False,
        "acetal": False,
        "tms": False,
    }

    deprotection_reactions = {
        "boc": False,
        "silyl": False,
        "benzyl": False,
        "acetal": False,
        "tms": False,
    }

    def dfs_traverse(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Boc protection/deprotection
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                protection_reactions["boc"] = True
                print("Found Boc protection reaction")

            if (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            ):
                deprotection_reactions["boc"] = True
                print("Found Boc deprotection reaction")

            # Silyl protection/deprotection
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                protection_reactions["silyl"] = True
                print("Found silyl protection reaction")

            if (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
                or checker.check_reaction("TMS deprotection from alkyne", rsmi)
            ):
                deprotection_reactions["silyl"] = True
                print("Found silyl deprotection reaction")

            # Benzyl protection/deprotection
            if (
                any(checker.check_fg("Primary alcohol", r) for r in reactants)
                and checker.check_fg("Ether", product)
                and "benzyl" in rsmi.lower()
            ):
                protection_reactions["benzyl"] = True
                print("Found benzyl protection reaction")

            if checker.check_reaction(
                "Hydroxyl benzyl deprotection", rsmi
            ) or checker.check_reaction("Carboxyl benzyl deprotection", rsmi):
                deprotection_reactions["benzyl"] = True
                print("Found benzyl deprotection reaction")

            # Acetal/ketal protection/deprotection
            if checker.check_reaction(
                "Aldehyde or ketone acetalization", rsmi
            ) or checker.check_reaction("Diol acetalization", rsmi):
                protection_reactions["acetal"] = True
                print("Found acetal/ketal protection reaction")

            if (
                checker.check_reaction("Acetal hydrolysis to diol", rsmi)
                or checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                or checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
            ):
                deprotection_reactions["acetal"] = True
                print("Found acetal/ketal deprotection reaction")

            # TMS protection/deprotection
            if checker.check_fg("TMS ether protective group", product) and any(
                checker.check_fg("Primary alcohol", r)
                or checker.check_fg("Secondary alcohol", r)
                or checker.check_fg("Tertiary alcohol", r)
                for r in reactants
            ):
                protection_reactions["tms"] = True
                print("Found TMS protection reaction")

            if any(checker.check_fg("TMS ether protective group", r) for r in reactants) and (
                checker.check_fg("Primary alcohol", product)
                or checker.check_fg("Secondary alcohol", product)
                or checker.check_fg("Tertiary alcohol", product)
            ):
                deprotection_reactions["tms"] = True
                print("Found TMS deprotection reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Count protection-deprotection cycles and different protecting groups used
    complete_cycles = 0
    different_groups = 0

    for group in protection_reactions:
        if protection_reactions[group] and deprotection_reactions[group]:
            complete_cycles += 1
        if protection_reactions[group] or deprotection_reactions[group]:
            different_groups += 1

    print(f"Complete protection-deprotection cycles: {complete_cycles}")
    print(f"Different protecting groups used: {different_groups}")

    # Return true if at least one complete protection-deprotection cycle
    # OR at least two different protecting groups are used
    return complete_cycles >= 1 or different_groups >= 2
