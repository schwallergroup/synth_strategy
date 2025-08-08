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
    This function detects a synthetic strategy involving late-stage amide formation
    (final step or penultimate step involves forming an amide bond).
    """
    # Track if the latest step is amide formation
    latest_step_is_amide_formation = False
    depth_of_amide_formation = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal latest_step_is_amide_formation, depth_of_amide_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for amide formation reactions
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Carboxylic acid to amide conversion",
                ]

                is_amide_formation = any(
                    checker.check_reaction(rxn_name, rsmi) for rxn_name in amide_formation_reactions
                )

                # If not detected by reaction type, check by functional groups
                if not is_amide_formation:
                    # Check if product has amide
                    product_has_amide = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )

                    # Check if reactants have acid derivatives and amines
                    has_acid_derivative = any(
                        checker.check_fg("Carboxylic acid", r)
                        or checker.check_fg("Ester", r)
                        or checker.check_fg("Acyl halide", r)
                        for r in reactants_smiles
                    )

                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants_smiles
                    )

                    is_amide_formation = has_acid_derivative and has_amine and product_has_amide

                if is_amide_formation:
                    print(f"Detected amide formation at depth {depth}")
                    if depth < depth_of_amide_formation:
                        depth_of_amide_formation = depth
                        # In retrosynthesis, depth 0 is target molecule, so reactions start at depth 1
                        # Consider it late-stage if it's within the first 2 steps (depth 1 or 2)
                        if depth <= 2:
                            latest_step_is_amide_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if latest_step_is_amide_formation:
        print("Detected late-stage amide formation strategy")
        return True
    else:
        print("No late-stage amide formation detected")
        return False
