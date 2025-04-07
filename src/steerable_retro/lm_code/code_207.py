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
    Detects if the synthesis route involves multiple ether formation reactions.
    """
    ether_formation_reactions = []

    def find_ether_formations(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]
            # Check if this is a Williamson ether synthesis or other ether formation
            if (
                checker.check_reaction("Williamson Ether Synthesis", rxn_smiles)
                or checker.check_reaction("Williamson Ether Synthesis (intra to epoxy)", rxn_smiles)
                or checker.check_reaction("Alcohol to ether", rxn_smiles)
                or checker.check_reaction("Mitsunobu aryl ether", rxn_smiles)
                or checker.check_reaction("Chan-Lam etherification", rxn_smiles)
                or checker.check_reaction("{Williamson ether}", rxn_smiles)
            ):
                ether_formation_reactions.append((rxn_smiles, depth))

            # Also check if product contains ether functional groups
            product = rxn_smiles.split(">")[-1]
            reactants = rxn_smiles.split(">")[0].split(".")

            # Check if ether is formed in the reaction (product has ether but reactants don't)
            if checker.check_fg("Ether", product):
                ether_in_reactants = any(checker.check_fg("Ether", r) for r in reactants)
                if not ether_in_reactants:
                    ether_formation_reactions.append((rxn_smiles, depth))

        for child in node.get("children", []):
            find_ether_formations(child, depth + 1)

    find_ether_formations(route)

    # Remove duplicates
    unique_reactions = set(rxn for rxn, _ in ether_formation_reactions)

    # Consider it a multiple ether formation strategy if at least 2 ether formations are found
    print(f"Ether formation reactions found: {len(unique_reactions)}")
    return len(unique_reactions) >= 2
