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
    This function detects a strategy where one hydroxyl of a diol is selectively
    functionalized (typically with a tosylate) while preserving the other hydroxyl.
    """
    # Initialize flags
    has_diol = False
    has_selective_functionalization = False

    def count_alcohols(smiles):
        """Helper function to count alcohol groups in a molecule"""
        count = 0
        if checker.check_fg("Primary alcohol", smiles):
            count += len(checker.get_fg_atom_indices("Primary alcohol", smiles))
        if checker.check_fg("Secondary alcohol", smiles):
            count += len(checker.get_fg_atom_indices("Secondary alcohol", smiles))
        if checker.check_fg("Tertiary alcohol", smiles):
            count += len(checker.get_fg_atom_indices("Tertiary alcohol", smiles))
        if checker.check_fg("Aromatic alcohol", smiles):
            count += len(checker.get_fg_atom_indices("Aromatic alcohol", smiles))
        return count

    def dfs_traverse(node, depth=0):
        nonlocal has_diol, has_selective_functionalization

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for diol in reactants
                for reactant in reactants:
                    if not reactant.strip():
                        continue

                    try:
                        # Count alcohol groups in reactant
                        alcohol_count_reactant = count_alcohols(reactant)

                        if alcohol_count_reactant >= 2:
                            has_diol = True
                            print(f"Detected diol at depth {depth}: {reactant}")

                            # Count alcohol groups in product
                            alcohol_count_product = count_alcohols(product)

                            # Check if exactly one alcohol was functionalized
                            if alcohol_count_product == alcohol_count_reactant - 1:
                                # Check for various functionalization reactions
                                functionalization_reactions = [
                                    "Formation of Sulfonic Esters",
                                    "Formation of Sulfonic Esters on TMS protected alcohol",
                                    "Alcohol to chloride_sulfonyl chloride",
                                    "Alcohol to chloride_SOCl2",
                                    "Alcohol to chloride_HCl",
                                    "Alcohol to chloride_Other",
                                    "Alcohol to triflate conversion",
                                    "Alcohol to bromide",
                                    "Alkyl iodides from alcohols",
                                    "Alkyl chlorides from alcohols",
                                    "Alkyl bromides from alcohols",
                                    "PBr3 and alcohol to alkyl bromide",
                                    "Appel reaction",
                                ]

                                for reaction_type in functionalization_reactions:
                                    if checker.check_reaction(reaction_type, rsmi):
                                        print(f"Detected reaction: {reaction_type}")
                                        has_selective_functionalization = True
                                        break

                                # If no specific reaction was detected, check for new leaving groups in product
                                if not has_selective_functionalization:
                                    leaving_groups = [
                                        "Tosylate",
                                        "Mesylate",
                                        "Triflate",
                                        "Primary halide",
                                        "Secondary halide",
                                        "Tertiary halide",
                                        "Aromatic halide",
                                    ]

                                    for lg in leaving_groups:
                                        if checker.check_fg(lg, product) and not checker.check_fg(
                                            lg, reactant
                                        ):
                                            print(f"Detected new leaving group: {lg}")
                                            has_selective_functionalization = True
                                            break

                                if has_selective_functionalization:
                                    print(f"Detected selective functionalization at depth {depth}")
                                    print(f"Reactant: {reactant}")
                                    print(f"Product: {product}")
                                    print(
                                        f"Alcohols before: {alcohol_count_reactant}, after: {alcohol_count_product}"
                                    )
                    except Exception as e:
                        print(f"Error processing molecule: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Strategy is present if both conditions are met
    strategy_present = has_diol and has_selective_functionalization

    if strategy_present:
        print("Detected selective diol functionalization strategy")
    else:
        print("Did not detect selective diol functionalization strategy")
        print(f"Has diol: {has_diol}")
        print(f"Has selective functionalization: {has_selective_functionalization}")

    return strategy_present
