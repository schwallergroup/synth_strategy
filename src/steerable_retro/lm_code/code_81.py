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
    Detects a synthesis route that includes an early-stage sulfonamide formation
    by coupling an amine with a sulfonyl chloride.
    """
    has_sulfonamide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_sulfonamide_formation

        if node["type"] == "reaction" and depth >= 2:  # Early stage (depth 2 or higher)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a sulfonamide formation reaction
                is_primary_sulfonamide = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                )
                is_secondary_sulfonamide = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )

                if is_primary_sulfonamide or is_secondary_sulfonamide:
                    print(f"Detected early-stage sulfonamide formation at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    has_sulfonamide_formation = True
                else:
                    # Fallback check using functional groups
                    has_amine = False
                    has_sulfonyl_chloride = False
                    has_sulfonamide_product = checker.check_fg("Sulfonamide", product)

                    for reactant in reactants:
                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            has_amine = True
                            print(f"Found amine reactant: {reactant}")
                        if checker.check_fg("Sulfonyl halide", reactant):
                            has_sulfonyl_chloride = True
                            print(f"Found sulfonyl chloride reactant: {reactant}")

                    if has_amine and has_sulfonyl_chloride and has_sulfonamide_product:
                        print(
                            f"Detected early-stage sulfonamide formation at depth {depth} using FG checks"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        has_sulfonamide_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_sulfonamide_formation
