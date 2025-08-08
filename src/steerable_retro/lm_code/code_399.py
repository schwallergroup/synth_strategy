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
    This function detects if the synthesis route uses late-stage amide formation
    as the final step in the synthesis.
    """
    amide_formation_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_at_late_stage

        print(f"Traversing node of type {node['type']} at depth {depth}")

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a late-stage reaction (depth 1 is the first reaction in retrosynthesis)
                if depth == 1:
                    print(f"Checking late-stage reaction: {rsmi}")

                    # Check for amide formation reactions using the checker
                    amide_formation_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Carboxylic acid with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Schotten-Baumann_amide",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                    ]

                    for reaction_name in amide_formation_reactions:
                        if checker.check_reaction(reaction_name, rsmi):
                            print(f"Detected late-stage amide formation: {reaction_name}")
                            amide_formation_at_late_stage = True
                            return

                    # If reaction checker didn't identify it, check for functional group changes
                    # Check if product has amide
                    if (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    ):
                        print(f"Product contains amide group: {product}")

                        # Check if reactants have necessary functional groups
                        has_carboxylic_acid = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants
                        )
                        has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                        has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                        has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants)

                        has_primary_amine = any(
                            checker.check_fg("Primary amine", r) for r in reactants
                        )
                        has_secondary_amine = any(
                            checker.check_fg("Secondary amine", r) for r in reactants
                        )
                        has_tertiary_amine = any(
                            checker.check_fg("Tertiary amine", r) for r in reactants
                        )
                        has_ammonia = (
                            any(checker.check_fg("Ammonia", r) for r in reactants)
                            if "Ammonia" in dir(checker)
                            else False
                        )

                        if (
                            has_carboxylic_acid or has_acyl_halide or has_ester or has_anhydride
                        ) and (has_primary_amine or has_secondary_amine or has_ammonia):
                            print(
                                "Detected late-stage amide formation based on functional group changes"
                            )
                            amide_formation_at_late_stage = True
                            return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)  # Increment depth as we go deeper in retrosynthesis

    # Start traversal from the root
    dfs_traverse(route)

    return amide_formation_at_late_stage
