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
    This function detects early-stage thiazole ring formation (depth ≥ 5)
    from thiourea and chloroacetyl compounds.
    """
    thiazole_formation_detected = False
    min_depth_for_early_stage = 5

    def dfs_traverse(node, depth=0):
        nonlocal thiazole_formation_detected

        if node["type"] == "reaction" and depth >= min_depth_for_early_stage:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a thiazole formation reaction
                is_thiazole_reaction = checker.check_reaction("thiazole", rsmi)

                # Check for thiourea in reactants
                thiourea_present = any(
                    checker.check_fg("Thiourea", reactant) for reactant in reactants
                )

                # Check for chloroacetyl compound in reactants (primary halide with adjacent carbonyl)
                chloroacetyl_present = any(
                    checker.check_fg("Primary halide", reactant)
                    and checker.check_fg("Ketone", reactant)
                    for reactant in reactants
                )

                # Check for thiazole ring in product
                thiazole_present = checker.check_ring("thiazole", product)

                print(f"Depth: {depth}, Thiazole reaction: {is_thiazole_reaction}")
                print(f"Thiourea present: {thiourea_present}")
                print(f"Chloroacetyl present: {chloroacetyl_present}")
                print(f"Thiazole present in product: {thiazole_present}")

                # Either detect the specific reaction type or check for the reactants and product
                if is_thiazole_reaction or (
                    thiourea_present and chloroacetyl_present and thiazole_present
                ):
                    print(f"Early-stage thiazole formation detected at depth {depth}")
                    thiazole_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return thiazole_formation_detected
