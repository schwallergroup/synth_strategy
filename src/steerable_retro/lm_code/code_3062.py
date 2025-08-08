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
    """Check if the route involves late-stage formation of an alcohol."""
    alcohol_formation_depths = []

    def dfs(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check if the product contains an alcohol group
            if (
                checker.check_fg("Primary alcohol", product)
                or checker.check_fg("Secondary alcohol", product)
                or checker.check_fg("Tertiary alcohol", product)
            ):

                # Check reactants to confirm this is an alcohol formation
                reactants = rsmi.split(">")[0].split(".")
                reactant_has_alcohol = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants
                )

                if not reactant_has_alcohol:
                    alcohol_formation_depths.append(depth)
                    print(f"Found alcohol formation at depth {depth}: {rsmi}")

                # Also check specific alcohol formation reactions
                elif (
                    checker.check_reaction("Reduction of aldehydes and ketones to alcohols", rsmi)
                    or checker.check_reaction("Grignard from aldehyde to alcohol", rsmi)
                    or checker.check_reaction("Grignard from ketone to alcohol", rsmi)
                    or checker.check_reaction("Reduction of ester to primary alcohol", rsmi)
                    or checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                    or checker.check_reaction(
                        "Reduction of carboxylic acid to primary alcohol", rsmi
                    )
                ):
                    alcohol_formation_depths.append(depth)
                    print(f"Found alcohol formation reaction at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)

    # Consider it late-stage if alcohol formation happens at a low depth (0-3)
    return any(depth <= 3 for depth in alcohol_formation_depths)
