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
    This function detects if C-N bond formation occurs in the last step (depth 0).
    """
    # List of C-N bond forming reactions
    cn_bond_forming_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Aminolysis of esters",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "Alkylation of amines",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Ugi reaction",
        "Goldberg coupling",
        "Ullmann-Goldberg Substitution amine",
    ]

    # Track if we found a late-stage C-N bond formation
    late_stage_cn_found = False

    def check_cn_bond_formation(reaction_node):
        """Check if this reaction forms a C-N bond"""
        if "metadata" not in reaction_node or "rsmi" not in reaction_node["metadata"]:
            return False

        rsmi = reaction_node["metadata"]["rsmi"]

        # Check if this is a known C-N bond forming reaction
        for reaction_type in cn_bond_forming_reactions:
            if checker.check_reaction(reaction_type, rsmi):
                print(f"Detected C-N bond forming reaction: {reaction_type}")
                return True

        # If not a known reaction type, check for C-N bond formation manually
        try:
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if product has nitrogen-containing functional groups
            product_has_n_fg = (
                checker.check_fg("Primary amine", product_part)
                or checker.check_fg("Secondary amine", product_part)
                or checker.check_fg("Tertiary amine", product_part)
                or checker.check_fg("Primary amide", product_part)
                or checker.check_fg("Secondary amide", product_part)
                or checker.check_fg("Tertiary amide", product_part)
                or checker.check_fg("Aniline", product_part)
            )

            if product_has_n_fg:
                # Split reactants and check if all have the same N-containing functional groups
                reactants = reactants_part.split(".")

                # Check if at least one reactant doesn't have N-containing functional groups
                # This would indicate N was introduced in this step
                for reactant in reactants:
                    if not (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Tertiary amine", reactant)
                        or checker.check_fg("Primary amide", reactant)
                        or checker.check_fg("Secondary amide", reactant)
                        or checker.check_fg("Tertiary amide", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        print("Detected C-N bond formation through functional group analysis")
                        return True
        except Exception as e:
            print(f"Error analyzing reaction: {e}")

        return False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cn_found

        print(f"Traversing node type: {node['type']} at depth {depth}")

        if node["type"] == "reaction":
            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                if check_cn_bond_formation(node):
                    print(f"Found late-stage C-N bond formation at depth {depth}")
                    late_stage_cn_found = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_cn_found
