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


def main(route, depth=0, sequence_count=0):
    """
    Detects if the synthesis route includes a sequence of aromatic substitution reactions.
    """
    if route["type"] == "reaction":
        try:
            rsmi = route["metadata"]["rsmi"]

            # Check if this is an aromatic substitution reaction
            is_aromatic_sub = False

            # Check for common aromatic substitution reactions
            if (
                checker.check_reaction("Aromatic chlorination", rsmi)
                or checker.check_reaction("Aromatic bromination", rsmi)
                or checker.check_reaction("Aromatic fluorination", rsmi)
                or checker.check_reaction("Aromatic iodination", rsmi)
                or checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                or checker.check_reaction("Friedel-Crafts alkylation", rsmi)
                or checker.check_reaction("Friedel-Crafts acylation", rsmi)
                or checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi)
            ):
                is_aromatic_sub = True

            if is_aromatic_sub:
                sequence_count += 1
                if sequence_count >= 2:  # At least 2 aromatic substitutions in sequence
                    print(
                        f"Aromatic substitution sequence detected with {sequence_count} reactions"
                    )
                    return True
        except Exception as e:
            print(f"Error analyzing aromatic substitution: {e}")

    # Recursively check children with updated sequence count
    for child in route.get("children", []):
        if aromatic_substitution_sequence_strategy(child, depth + 1, sequence_count):
            return True

    return False
