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
    Detects a synthesis where the final step is a deprotection reaction,
    specifically focusing on Boc deprotection.
    """
    has_late_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_deprotection

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a Boc deprotection reaction
                if checker.check_reaction("Boc amine deprotection", rsmi):
                    print(f"Found Boc amine deprotection reaction at depth {depth}: {rsmi}")
                    has_late_deprotection = True
                    return

                # Check for other common deprotection reactions
                deprotection_reactions = [
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine",
                    "Hydroxyl benzyl deprotection",
                    "Carboxyl benzyl deprotection",
                    "Cleavage of methoxy ethers to alcohols",
                    "Cleavage of alkoxy ethers to alcohols",
                    "COOH ethyl deprotection",
                    "N-glutarimide deprotection",
                    "Phthalimide deprotection",
                    "TMS deprotection from alkyne",
                    "Alcohol deprotection from silyl ethers",
                    "Alcohol deprotection from silyl ethers (double)",
                    "Alcohol deprotection from silyl ethers (diol)",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                ]

                for reaction_type in deprotection_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found {reaction_type} reaction at depth {depth}: {rsmi}")
                        has_late_deprotection = True
                        return

                # Manual check for common protecting groups
                # Boc group
                if any(
                    checker.check_fg("Boc", r) for r in reactants_smiles
                ) and not checker.check_fg("Boc", product_smiles):
                    print(f"Found Boc group removal at depth {depth} (manual check): {rsmi}")
                    has_late_deprotection = True
                    return

                # Silyl protecting groups
                if any(
                    checker.check_fg("TMS ether protective group", r)
                    or checker.check_fg("Silyl protective group", r)
                    for r in reactants_smiles
                ) and not (
                    checker.check_fg("TMS ether protective group", product_smiles)
                    or checker.check_fg("Silyl protective group", product_smiles)
                ):
                    print(f"Found silyl group removal at depth {depth} (manual check): {rsmi}")
                    has_late_deprotection = True
                    return

                # Acetal/Ketal deprotection
                if any(
                    checker.check_fg("Acetal/Ketal", r) for r in reactants_smiles
                ) and not checker.check_fg("Acetal/Ketal", product_smiles):
                    if checker.check_fg("Aldehyde", product_smiles) or checker.check_fg(
                        "Ketone", product_smiles
                    ):
                        print(
                            f"Found acetal/ketal deprotection at depth {depth} (manual check): {rsmi}"
                        )
                        has_late_deprotection = True
                        return

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_late_deprotection
