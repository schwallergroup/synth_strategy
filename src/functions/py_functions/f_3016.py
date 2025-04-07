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
    This function detects multiple C-O bond formations in the synthesis.
    Looks for at least 2 reactions forming C-O bonds.
    """
    co_bond_formation_count = 0

    def dfs_traverse(node):
        nonlocal co_bond_formation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for C-O bond forming reactions using the checker function
                co_forming_reactions = [
                    "Williamson Ether Synthesis",
                    "Williamson Ether Synthesis (intra to epoxy)",
                    "Esterification of Carboxylic Acids",
                    "Formation of Sulfonic Esters",
                    "Formation of Sulfonic Esters on TMS protected alcohol",
                    "Alcohol protection with silyl ethers",
                    "Oxidative esterification of primary alcohols",
                    "Transesterification",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Schotten-Baumann to ester",
                    "Mitsunobu aryl ether",
                    "Mitsunobu esterification",
                    "Mitsunobu aryl ether (intramolecular)",
                    "Chan-Lam alcohol",
                    "Chan-Lam etherification",
                    "Acetic anhydride and alcohol to ester",
                ]

                for reaction_type in co_forming_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        co_bond_formation_count += 1
                        print(f"Detected C-O bond formation: {reaction_type}")
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Current count: {co_bond_formation_count}")
                        break

                # If no specific reaction type matched, check for general C-O bond formation
                if not any(
                    checker.check_reaction(r, rsmi) for r in co_forming_reactions
                ):
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for alcohol/phenol in reactants
                    has_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Phenol", r)
                        for r in reactants
                    )

                    # Check for ether/ester in product
                    has_ether_ester = checker.check_fg(
                        "Ether", product
                    ) or checker.check_fg("Ester", product)

                    if has_alcohol and has_ether_ester:
                        co_bond_formation_count += 1
                        print(
                            f"Detected general C-O bond formation (alcohol → ether/ester)"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Current count: {co_bond_formation_count}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total C-O bond formations found: {co_bond_formation_count}")
    return co_bond_formation_count >= 2  # Return True if at least 2 C-O bond formations
