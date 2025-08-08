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
    This function detects a synthetic strategy featuring amide formation from an acyl chloride
    and amine in the early stage of the synthesis.
    """
    has_amide_formation = False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_amide_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Try to get depth from metadata if available
            depth = current_depth
            if "metadata" in node and "ID" in node["metadata"]:
                try:
                    id_string = node["metadata"]["ID"]
                    if "Depth: " in id_string:
                        depth_str = id_string.split("Depth: ")[1].split()[0]
                        depth = int(depth_str)
                except Exception as e:
                    print(f"Error extracting depth: {e}")

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Early stage is considered depth >= 2
            if depth >= 2:
                print(f"Checking early-stage reaction at depth {depth}: {rsmi}")

                # Check if this is an amide formation from acyl chloride using reaction checkers
                is_amide_formation = (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                )

                if is_amide_formation:
                    print(f"Reaction checker identified amide formation: {rsmi}")

                # Additional verification by checking functional groups
                has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                has_amine = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    for r in reactants
                )
                has_amide = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                print(
                    f"FG checks - Acyl halide: {has_acyl_halide}, Amine: {has_amine}, Amide: {has_amide}"
                )

                if is_amide_formation and has_acyl_halide and has_amine and has_amide:
                    print("Detected amide formation from acyl chloride in early step")
                    has_amide_formation = True
                elif has_acyl_halide and has_amine and has_amide:
                    # Even if reaction checker failed, if we have the right FGs, it's likely an amide formation
                    print("Detected likely amide formation based on functional groups")
                    has_amide_formation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_amide_formation
