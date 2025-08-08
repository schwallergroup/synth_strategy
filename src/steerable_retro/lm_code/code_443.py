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
    This function detects a synthetic strategy involving amine to nitrile conversion
    in the early stages of synthesis.
    """
    amine_to_nitrile_found = False
    # Define what constitutes "early stages" of synthesis
    EARLY_STAGE_THRESHOLD = 3

    def dfs_traverse(node, depth=0):
        nonlocal amine_to_nitrile_found

        if node["type"] == "reaction" and depth < EARLY_STAGE_THRESHOLD:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amine in reactants and nitrile in product (forward direction)
                amine_in_reactants = any(
                    checker.check_fg("Primary amine", reactant) for reactant in reactants
                )
                nitrile_in_product = checker.check_fg("Nitrile", product)

                # Check for nitrile in reactants and amine in product (reverse direction)
                nitrile_in_reactants = any(
                    checker.check_fg("Nitrile", reactant) for reactant in reactants
                )
                amine_in_product = checker.check_fg("Primary amine", product)

                # Forward direction: amine to nitrile
                if amine_in_reactants and nitrile_in_product:
                    print(f"Found potential amine to nitrile conversion at depth {depth}")
                    # Check for common reactions or intermediates in this pathway
                    if any(
                        checker.check_reaction(rxn_type, rsmi)
                        for rxn_type in ["Amine to azide", "Azide to nitrile", "Dehydration"]
                    ):
                        print(f"Confirmed amine to nitrile conversion reaction at depth {depth}")
                        amine_to_nitrile_found = True
                    else:
                        # If no specific reaction is identified but pattern exists, accept it
                        print(f"Direct amine to nitrile pattern found at depth {depth}")
                        amine_to_nitrile_found = True

                # Reverse direction: nitrile to amine (in retrosynthesis means amine to nitrile in forward)
                if nitrile_in_reactants and amine_in_product:
                    print(f"Found potential nitrile to amine conversion at depth {depth}")
                    if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                        print(f"Confirmed nitrile to amine reduction reaction at depth {depth}")
                        amine_to_nitrile_found = True
                    else:
                        # If no specific reaction is identified but pattern exists, accept it
                        print(f"Direct nitrile to amine pattern found at depth {depth}")
                        amine_to_nitrile_found = True

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return amine_to_nitrile_found
