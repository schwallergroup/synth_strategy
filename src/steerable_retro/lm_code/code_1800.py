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
    This function detects if the synthetic route involves a late-stage reductive amination.
    Late-stage is defined as occurring at depth 0 or 1 in the synthesis tree.
    """
    has_late_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_reductive_amination

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for reductive amination reaction types
                if (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                ):

                    print(f"Found late-stage reductive amination at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    has_late_reductive_amination = True
                    return  # Exit early once found

                # Fallback check if reaction type checker fails
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde or ketone in reactants
                has_carbonyl = False
                for r in reactants:
                    if (
                        checker.check_fg("Aldehyde", r)
                        or checker.check_fg("Ketone", r)
                        or checker.check_fg("Formaldehyde", r)
                    ):
                        has_carbonyl = True
                        break

                # Check for amine in reactants
                has_amine = False
                for r in reactants:
                    if checker.check_fg("Primary amine", r) or checker.check_fg(
                        "Secondary amine", r
                    ):
                        has_amine = True
                        break

                # Check for new amine in product that wasn't in reactants
                if has_carbonyl and has_amine:
                    # Additional check to confirm this is likely a reductive amination
                    has_carbonyl_in_product = (
                        checker.check_fg("Aldehyde", product)
                        or checker.check_fg("Ketone", product)
                        or checker.check_fg("Formaldehyde", product)
                    )

                    if not has_carbonyl_in_product:
                        print(
                            f"Found potential late-stage reductive amination at depth {depth} (fallback detection)"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        has_late_reductive_amination = True
                        return  # Exit early once found

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_late_reductive_amination
