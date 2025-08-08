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
    This function detects a synthetic strategy involving carbonyl reduction
    (ketone or aldehyde to secondary alcohol).
    """
    has_carbonyl_reduction = False

    def dfs_traverse(node):
        nonlocal has_carbonyl_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Get depth from metadata
            depth = -1
            if "ID" in node["metadata"]:
                depth_str = node["metadata"]["ID"]
                if "Depth:" in depth_str:
                    try:
                        depth = int(depth_str.split("Depth:")[1].split()[0])
                    except:
                        pass

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check for carbonyl reduction using specific reaction types
            if checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi):
                print(f"Found specific ketone reduction reaction at depth {depth}")
                has_carbonyl_reduction = True
            elif checker.check_reaction("Reduction of aldehydes and ketones to alcohols", rsmi):
                # Verify it's specifically carbonyl to secondary alcohol
                carbonyl_in_reactants = any(
                    checker.check_fg("Ketone", r) or checker.check_fg("Aldehyde", r)
                    for r in reactants
                )
                sec_alcohol_in_product = checker.check_fg("Secondary alcohol", product)
                sec_alcohol_in_reactants = any(
                    checker.check_fg("Secondary alcohol", r) for r in reactants
                )

                if (
                    carbonyl_in_reactants
                    and sec_alcohol_in_product
                    and not sec_alcohol_in_reactants
                ):
                    print(
                        f"Found general reduction reaction with carbonyl to secondary alcohol at depth {depth}"
                    )
                    has_carbonyl_reduction = True
            else:
                # Alternative check: look for carbonyl in reactants and secondary alcohol in product
                reactants_str = ".".join(reactants)
                for reactant in reactants:
                    if (
                        checker.check_fg("Ketone", reactant)
                        or checker.check_fg("Aldehyde", reactant)
                    ) and not checker.check_fg("Secondary alcohol", reactant):
                        if checker.check_fg("Secondary alcohol", product) and not (
                            checker.check_fg("Ketone", product)
                            and checker.check_fg("Aldehyde", product)
                        ):
                            print(
                                f"Detected carbonyl reduction by functional group analysis at depth {depth}"
                            )
                            has_carbonyl_reduction = True
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Final result: carbonyl reduction detected = {has_carbonyl_reduction}")
    return has_carbonyl_reduction
