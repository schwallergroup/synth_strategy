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
    This function detects if the synthesis involves a late-stage amide coupling
    (in the final or penultimate step).
    """
    has_late_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_amide_coupling

        if node["type"] == "reaction" and depth <= 1:  # Only check final or penultimate steps
            try:
                rsmi = node["metadata"]["rsmi"]

                # Check for amide coupling reactions using the checker function
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Ester with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                ]

                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected late-stage amide coupling reaction: {reaction_type} at depth {depth}"
                        )
                        has_late_amide_coupling = True
                        break

                # If no specific reaction type matched, check for functional group patterns
                if not has_late_amide_coupling:
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    # Check if product contains amide group
                    has_amide_product = (
                        checker.check_fg("Primary amide", product_part)
                        or checker.check_fg("Secondary amide", product_part)
                        or checker.check_fg("Tertiary amide", product_part)
                    )

                    if has_amide_product:
                        # Check if reactants contain carboxylic acid and amine
                        reactants = reactants_part.split(".")
                        has_carboxylic_acid = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants
                        )
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Tertiary amine", r)
                            for r in reactants
                        )

                        if has_carboxylic_acid and has_amine:
                            print(
                                f"Detected late-stage amide coupling via functional groups at depth {depth}"
                            )
                            has_late_amide_coupling = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Late-stage amide coupling: {has_late_amide_coupling}")
    return has_late_amide_coupling
