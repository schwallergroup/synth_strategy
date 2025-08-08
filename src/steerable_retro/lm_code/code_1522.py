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
    Detects if azide functional group is introduced in the late stages of synthesis
    (within the first half of the synthesis depth).
    """
    max_depth = 0
    azide_introduction_depth = None

    # First pass to determine maximum depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to find where azide is introduced
    def find_azide_introduction(node, depth=0):
        nonlocal azide_introduction_depth

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check if product has azide but reactants don't
                if checker.check_fg("Azide", product_smiles):
                    if not any(
                        checker.check_fg("Azide", reactant) for reactant in reactants_smiles
                    ):
                        # In retrosynthesis, we're going backward, so this is where azide is introduced
                        if azide_introduction_depth is None or depth < azide_introduction_depth:
                            azide_introduction_depth = depth
                            print(f"Azide introduced at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            find_azide_introduction(child, depth + 1)

    # Run the traversals
    find_max_depth(route)
    find_azide_introduction(route)

    # Check if azide is introduced in the first half of synthesis
    if azide_introduction_depth is not None:
        # In retrosynthesis, late-stage means higher depth values (closer to starting materials)
        is_late_stage = azide_introduction_depth <= (max_depth / 2)
        print(f"Max depth: {max_depth}, Azide introduction depth: {azide_introduction_depth}")
        print(f"Azide introduction is {'late-stage' if is_late_stage else 'early-stage'}")
        return is_late_stage
    else:
        print("No azide introduction found in the synthesis route")

    return False
