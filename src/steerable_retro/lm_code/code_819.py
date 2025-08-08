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
    This function detects if carbamate formation/cleavage occurs in the late stage of synthesis.
    """
    late_stage_carbamate = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_carbamate

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
            and depth <= 2
        ):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for carbamate-related reactions directly
                carbamate_reactions = [
                    "Boc amine protection",
                    "Boc amine deprotection",
                    "Boc amine protection explicit",
                    "Boc amine protection with Boc anhydride",
                    "Boc amine protection (ethyl Boc)",
                    "Boc amine protection of secondary amine",
                    "Boc amine protection of primary amine",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine",
                    "Urea synthesis via isocyanate and primary amine",
                    "Urea synthesis via isocyanate and secondary amine",
                    "Urea synthesis via isocyanate and diazo",
                    "Urea synthesis via isocyanate and sulfonamide",
                ]

                for reaction_name in carbamate_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        print(
                            f"Detected carbamate reaction {reaction_name} at depth {depth}: {rsmi}"
                        )
                        late_stage_carbamate = True
                        return

                # Check for functional group changes if specific reaction check failed
                # Check for carbamate-related functional groups in reactants
                carbamate_fgs = ["Urea", "Carbonic Ester", "Carbonic Acid"]
                carbamate_in_reactants = any(
                    any(checker.check_fg(fg, r) for fg in carbamate_fgs) for r in reactants_smiles
                )
                carbamate_in_product = any(
                    checker.check_fg(fg, product_smiles) for fg in carbamate_fgs
                )

                # Check for amine functional groups in reactants and products
                amine_fgs = ["Primary amine", "Secondary amine", "Tertiary amine"]
                amine_in_reactants = any(
                    any(checker.check_fg(fg, r) for fg in amine_fgs) for r in reactants_smiles
                )
                amine_in_product = any(checker.check_fg(fg, product_smiles) for fg in amine_fgs)

                # Check for isocyanate which is often involved in carbamate formation
                isocyanate_in_reactants = any(
                    checker.check_fg("Isocyanate", r) for r in reactants_smiles
                )

                # Check for carbamate formation or cleavage
                if (
                    (amine_in_reactants and carbamate_in_product)
                    or (carbamate_in_reactants and amine_in_product)
                    or (isocyanate_in_reactants and (carbamate_in_product or amine_in_product))
                ):
                    print(f"Detected carbamate formation/cleavage at depth {depth}: {rsmi}")
                    late_stage_carbamate = True
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_carbamate
