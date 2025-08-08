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
    This function detects a synthetic strategy involving late-stage introduction
    of a trifluoromethyl group.
    """
    has_late_trifluoromethyl_introduction = False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_late_trifluoromethyl_introduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                # Extract depth - try from ID first, otherwise use current_depth
                depth = -1
                if "ID" in node["metadata"]:
                    try:
                        depth_str = node["metadata"]["ID"].split("Depth: ")[-1].split()[0]
                        depth = int(depth_str)
                    except (IndexError, ValueError):
                        depth = current_depth
                else:
                    depth = current_depth

                print(f"Examining reaction at depth {depth}")

                # Consider reactions at depth 0, 1, or 2 as late-stage
                if depth <= 2:
                    rsmi = node["metadata"]["rsmi"]
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check for trifluoromethyl group in product
                    product_has_cf3 = checker.check_fg("Trifluoro group", product_smiles)
                    print(f"Product has CF3: {product_has_cf3}")

                    if product_has_cf3:
                        # Check if any reactant doesn't have CF3 (indicating introduction)
                        reactants_with_cf3 = [
                            checker.check_fg("Trifluoro group", reactant)
                            for reactant in reactants_smiles
                        ]
                        print(f"Reactants with CF3: {reactants_with_cf3}")

                        # If at least one reactant doesn't have CF3, and product has CF3,
                        # then CF3 was introduced in this reaction
                        if not all(reactants_with_cf3) and any(reactants_with_cf3):
                            print(
                                f"Detected late-stage introduction of trifluoromethyl group at depth {depth}"
                            )
                            has_late_trifluoromethyl_introduction = True

                        # If no reactants have CF3 but product does, it's definitely introduced
                        if not any(reactants_with_cf3):
                            print(
                                f"Detected late-stage introduction of trifluoromethyl group at depth {depth}"
                            )
                            has_late_trifluoromethyl_introduction = True

                        # Check for specific CF3 introduction reactions
                        if checker.check_reaction("Fluorination", rsmi) or checker.check_reaction(
                            "Trifluoromethylation", rsmi
                        ):
                            print(f"Detected specific CF3 introduction reaction at depth {depth}")
                            has_late_trifluoromethyl_introduction = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {has_late_trifluoromethyl_introduction}")

    return has_late_trifluoromethyl_introduction
