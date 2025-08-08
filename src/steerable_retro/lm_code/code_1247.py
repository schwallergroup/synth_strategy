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
    Detects if a cycloalkyl group is incorporated at a late stage in the synthesis.
    Late stage is defined as being in the first 1/3 of the synthesis depth.
    """
    max_depth = 0
    cycloalkyl_incorporations = []

    def find_cycloalkyl_incorporations(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]
            reactants = rxn_smiles.split(">")[0].split(".")
            product = rxn_smiles.split(">")[-1]

            # Check if any reactant contains a cycloalkyl group that's transferred to the product
            has_cycloalkyl_transfer = False
            for reactant in reactants:
                if (
                    checker.check_ring("cyclopropane", reactant)
                    or checker.check_ring("cyclobutane", reactant)
                    or checker.check_ring("cyclopentane", reactant)
                    or checker.check_ring("cyclohexane", reactant)
                    or checker.check_ring("cycloheptane", reactant)
                    or checker.check_ring("cyclooctane", reactant)
                ):

                    # Verify the cycloalkyl group is in the product
                    if (
                        checker.check_ring("cyclopropane", product)
                        or checker.check_ring("cyclobutane", product)
                        or checker.check_ring("cyclopentane", product)
                        or checker.check_ring("cyclohexane", product)
                        or checker.check_ring("cycloheptane", product)
                        or checker.check_ring("cyclooctane", product)
                    ):
                        has_cycloalkyl_transfer = True
                        break

            if has_cycloalkyl_transfer:
                cycloalkyl_incorporations.append((rxn_smiles, depth))

        for child in node.get("children", []):
            find_cycloalkyl_incorporations(child, depth + 1)

    find_cycloalkyl_incorporations(route)

    if not cycloalkyl_incorporations or max_depth == 0:
        return False

    # Sort by depth (ascending)
    cycloalkyl_incorporations.sort(key=lambda x: x[1])

    # Check if any cycloalkyl incorporation happens in the first 1/3 of the synthesis
    late_stage_threshold = max_depth / 3
    for _, depth in cycloalkyl_incorporations:
        if depth <= late_stage_threshold:
            print(
                f"Late-stage cycloalkyl incorporation found at depth {depth} (threshold: {late_stage_threshold:.1f})"
            )
            return True

    print(
        f"No late-stage cycloalkyl incorporation found. Earliest at depth {cycloalkyl_incorporations[0][1]} (threshold: {late_stage_threshold:.1f})"
    )
    return False
