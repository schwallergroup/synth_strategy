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
    This function detects a synthetic strategy involving nitration followed by
    reduction of the nitro group to an amine.
    """
    # Track if we found the relevant reactions
    found_nitration = False
    found_reduction = False
    nitro_depth = None
    reduction_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_nitration, found_reduction, nitro_depth, reduction_depth

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Extract reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitration (introduction of nitro group)
                is_nitration = (
                    checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                    or checker.check_reaction("Non-aromatic nitration with HNO3", rsmi)
                )

                # Additional check: product should have nitro group that wasn't in reactants
                if (
                    is_nitration
                    or any(not checker.check_fg("Nitro group", r) for r in reactants)
                    and checker.check_fg("Nitro group", product)
                ):
                    found_nitration = True
                    nitro_depth = depth
                    print(f"Found nitration at depth {depth}")

                # Check for nitro reduction
                is_reduction = checker.check_reaction("Reduction of nitro groups to amines", rsmi)

                # Additional check: reactant should have nitro group and product should have amine
                if is_reduction or (
                    any(checker.check_fg("Nitro group", r) for r in reactants)
                    and not checker.check_fg("Nitro group", product)
                    and checker.check_fg("Aniline", product)
                    or checker.check_fg("Primary amine", product)
                ):
                    found_reduction = True
                    reduction_depth = depth
                    print(f"Found nitro reduction at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitration found: {found_nitration}, depth: {nitro_depth}")
    print(f"Reduction found: {found_reduction}, depth: {reduction_depth}")

    # Return True if both reactions were found in the correct sequence
    # (reduction should be at a lower depth than nitration in the tree)
    return (
        found_nitration
        and found_reduction
        and (
            reduction_depth is not None
            and nitro_depth is not None
            and reduction_depth < nitro_depth
        )
    )
