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
    This function detects a synthesis strategy involving late-stage nitro reduction.
    """
    has_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro reduction reaction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Found nitro reduction reaction at depth {depth}: {rsmi}")

                # Verify reactants have nitro groups and product has primary amines
                nitro_in_reactants = any(
                    checker.check_fg("Nitro group", reactant) for reactant in reactants
                )
                amine_in_product = checker.check_fg("Primary amine", product)

                if nitro_in_reactants and amine_in_product:
                    # Check if this is a late-stage reaction (depth 0, 1, or 2)
                    if depth <= 2:
                        print(f"Confirmed late-stage nitro reduction at depth {depth}")
                        has_nitro_reduction = True
                else:
                    print(
                        f"Reaction is nitro reduction but functional groups don't match: nitro={nitro_in_reactants}, amine={amine_in_product}"
                    )

            # Alternative check in case the reaction checker fails
            elif not has_nitro_reduction:
                for reactant in reactants:
                    if checker.check_fg("Nitro group", reactant) and checker.check_fg(
                        "Primary amine", product
                    ):
                        print(
                            f"Potential nitro reduction detected at depth {depth} (not confirmed by reaction checker)"
                        )
                        # Only set to true if late-stage and we're confident it's a nitro reduction
                        if depth <= 2 and not any(
                            checker.check_fg("Primary amine", r) for r in reactants
                        ):
                            print(
                                f"Setting as late-stage nitro reduction based on functional group analysis"
                            )
                            has_nitro_reduction = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: has_nitro_reduction = {has_nitro_reduction}")
    return has_nitro_reduction
