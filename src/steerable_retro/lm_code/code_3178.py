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
    This function detects amide bond formation from carboxylic acid and amine.
    """
    amide_formation_detected = False

    def dfs_traverse(node):
        nonlocal amide_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an amide formation reaction using the checker function
            amide_formation_reactions = [
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Carboxylic acid with primary amine to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Ester with ammonia to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Acyl chloride with ammonia to amide",
            ]

            is_amide_formation = any(
                checker.check_reaction(reaction_name, rsmi)
                for reaction_name in amide_formation_reactions
            )

            if is_amide_formation:
                print(f"Amide formation reaction detected: {rsmi}")
                amide_formation_detected = True
            else:
                # Fallback method: check for functional group transformation
                has_carboxylic_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                has_secondary_amine = any(checker.check_fg("Secondary amine", r) for r in reactants)
                has_amide = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                if has_carboxylic_acid and (has_primary_amine or has_secondary_amine) and has_amide:
                    print(f"Amide formation detected through functional group analysis: {rsmi}")
                    print(f"  - Carboxylic acid in reactants: {has_carboxylic_acid}")
                    print(f"  - Primary amine in reactants: {has_primary_amine}")
                    print(f"  - Secondary amine in reactants: {has_secondary_amine}")
                    print(f"  - Amide in product: {has_amide}")
                    amide_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: Amide formation detected = {amide_formation_detected}")
    return amide_formation_detected
