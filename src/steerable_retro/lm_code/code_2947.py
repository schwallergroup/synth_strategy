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
    Detects a strategy involving sequential formation of multiple heterocycles,
    specifically pyrazole formation followed by thiazole formation.
    """
    pyrazole_formation = False
    thiazole_formation = False
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal pyrazole_formation, thiazole_formation, reaction_sequence

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for pyrazole formation
                if checker.check_reaction("pyrazole", rsmi) or (
                    checker.check_ring("pyrazole", product)
                    and not any(checker.check_ring("pyrazole", r) for r in reactants)
                ):
                    pyrazole_formation = True
                    reaction_sequence.append(("pyrazole_formation", depth))
                    print(f"Detected pyrazole formation at depth {depth}, RSMI: {rsmi}")

                # Check for thiazole formation
                if checker.check_reaction("thiazole", rsmi) or (
                    checker.check_ring("thiazole", product)
                    and not any(checker.check_ring("thiazole", r) for r in reactants)
                ):
                    thiazole_formation = True
                    reaction_sequence.append(("thiazole_formation", depth))
                    print(f"Detected thiazole formation at depth {depth}, RSMI: {rsmi}")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Pyrazole formation detected: {pyrazole_formation}")
    print(f"Thiazole formation detected: {thiazole_formation}")
    print(f"Reaction sequence: {reaction_sequence}")

    # Check if both heterocycles were formed
    if pyrazole_formation and thiazole_formation:
        # Find depths of each reaction
        pyrazole_depth = next(
            (depth for rxn, depth in reaction_sequence if rxn == "pyrazole_formation"), -1
        )
        thiazole_depth = next(
            (depth for rxn, depth in reaction_sequence if rxn == "thiazole_formation"), -1
        )

        # Higher depth = earlier in synthesis (retrosynthetic perspective)
        correct_sequence = pyrazole_depth > thiazole_depth

        if correct_sequence:
            print(
                "Detected sequential heterocycle formation strategy: pyrazole followed by thiazole"
            )
            return True
        else:
            print(
                f"Heterocycles formed in wrong order: pyrazole at depth {pyrazole_depth}, thiazole at depth {thiazole_depth}"
            )

    return False
