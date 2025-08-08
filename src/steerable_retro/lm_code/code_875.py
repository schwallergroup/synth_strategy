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
    This function detects the strategy of sequential interconversion of leaving groups
    (e.g., azide to chloride or other leaving group transformations).
    """
    leaving_group_changes = 0

    # List of leaving groups to check
    leaving_groups = [
        "Azide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Alkenyl halide",
        "Aromatic halide",
        "Triflate",
        "Mesylate",
        "Tosylate",
    ]

    # List of specific reactions that involve leaving group interconversions
    lg_conversion_reactions = [
        "Formation of Azides from halogens",
        "Formation of Azides from boronic acids",
        "Alcohol to chloride_sulfonyl chloride",
        "Alcohol to chloride_SOCl2",
        "Alcohol to chloride_CHCl3",
        "Alcohol to chloride_CH2Cl2",
        "Alcohol to chloride_PCl5_ortho",
        "Alcohol to chloride_POCl3_ortho",
        "Alcohol to chloride_POCl3_para",
        "Alcohol to chloride_POCl3",
        "Alcohol to chloride_HCl",
        "Alcohol to chloride_Salt",
        "Alcohol to chloride_Other",
        "Alcohol to triflate conversion",
        "Primary amine to fluoride",
        "Primary amine to chloride",
        "Primary amine to bromide",
        "Primary amine to iodide",
        "Amine to azide",
        "Appel reaction",
        "PBr3 and alcohol to alkyl bromide",
        "Finkelstein reaction",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal leaving_group_changes

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a known leaving group interconversion reaction
            for reaction_type in lg_conversion_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected leaving group interconversion reaction: {reaction_type}")
                    leaving_group_changes += 1
                    break
            else:
                # If not a known reaction type, check for changes in leaving groups
                reactant_groups = {}
                product_groups = {}

                # Check reactants for leaving groups
                for reactant in reactants:
                    for lg in leaving_groups:
                        if checker.check_fg(lg, reactant):
                            if lg not in reactant_groups:
                                reactant_groups[lg] = 0
                            reactant_groups[lg] += 1

                # Check product for leaving groups
                for lg in leaving_groups:
                    if checker.check_fg(lg, product):
                        if lg not in product_groups:
                            product_groups[lg] = 0
                        product_groups[lg] += 1

                # If there's a change in leaving groups, count it
                if reactant_groups and product_groups and reactant_groups != product_groups:
                    # Additional check: ensure we're not just adding a leaving group
                    # but actually replacing one with another
                    if any(
                        lg in reactant_groups and lg not in product_groups for lg in leaving_groups
                    ) and any(
                        lg not in reactant_groups and lg in product_groups for lg in leaving_groups
                    ):
                        print(
                            f"Detected leaving group change: {reactant_groups} -> {product_groups}"
                        )
                        leaving_group_changes += 1

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    result = leaving_group_changes >= 1
    print(
        f"Leaving group interconversion strategy detected: {result} (changes: {leaving_group_changes})"
    )
    return result
