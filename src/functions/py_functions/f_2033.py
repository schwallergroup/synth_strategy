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
    Detects if the synthesis includes oxidation of a terminal alkene to an aldehyde in the middle of the synthesis.
    In retrosynthesis, this appears as a transformation from aldehyde to alkene (e.g., Wittig reaction).
    """
    found_alkene_oxidation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_alkene_oxidation

        if node["type"] == "reaction" and 1 <= depth <= 4:  # Mid-synthesis reactions
            print(f"Examining reaction at depth {depth}")
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                print(f"Reaction SMILES: {rsmi}")

                # In retrosynthesis, alkene oxidation appears as aldehyde to alkene transformation
                # Check for Wittig reactions which are the retrosynthetic equivalent of alkene oxidation
                is_wittig = (
                    checker.check_reaction(
                        "Wittig reaction with triphenylphosphorane", rsmi
                    )
                    or checker.check_reaction("Wittig with Phosphonium", rsmi)
                    or checker.check_reaction("Wittig", rsmi)
                )
                print(f"Is Wittig reaction: {is_wittig}")

                # Check for Julia Olefination which is another aldehyde to alkene transformation
                is_julia = checker.check_reaction("Julia Olefination", rsmi)
                print(f"Is Julia Olefination: {is_julia}")

                # Check functional groups: in retrosynthesis, we need aldehyde in reactants, alkene in product
                has_aldehyde = False
                for reactant in reactants:
                    if checker.check_fg("Aldehyde", reactant) or checker.check_fg(
                        "Formaldehyde", reactant
                    ):
                        print(f"Found aldehyde in reactant: {reactant}")
                        has_aldehyde = True
                        break

                has_alkene = (
                    checker.check_fg("Vinyl", product)
                    or checker.check_fg("Ethylene", product)
                    or checker.check_fg("Alkene", product)
                    or checker.check_fg("Allene", product)
                    or checker.check_fg("Conjugated diene", product)
                )
                print(f"Has alkene in product: {has_alkene}")

                # If we find a Wittig/Julia reaction or the pattern of aldehyde → alkene
                if is_wittig or is_julia or (has_aldehyde and has_alkene):
                    print(f"Found alkene oxidation strategy at depth {depth}")
                    found_alkene_oxidation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {found_alkene_oxidation}")
    return found_alkene_oxidation
