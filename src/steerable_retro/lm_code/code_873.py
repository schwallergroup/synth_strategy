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
    This function detects the strategy of converting amines to azides and potentially back,
    which is a common approach for amine functionalization or protection.
    """
    amine_to_azide_conversion = False
    azide_to_amine_conversion = False

    def dfs_traverse(node):
        nonlocal amine_to_azide_conversion, azide_to_amine_conversion

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amine to azide conversion
            reactants_have_amine = any(
                checker.check_fg("Primary amine", smi)
                or checker.check_fg("Secondary amine", smi)
                or checker.check_fg("Tertiary amine", smi)
                for smi in reactants
            )

            product_has_azide = checker.check_fg("Azide", product)

            # Check for azide to amine conversion
            reactants_have_azide = any(checker.check_fg("Azide", smi) for smi in reactants)

            product_has_amine = (
                checker.check_fg("Primary amine", product)
                or checker.check_fg("Secondary amine", product)
                or checker.check_fg("Tertiary amine", product)
            )

            # Check for specific reactions
            amine_to_azide_reaction = checker.check_reaction("Amine to azide", rsmi)
            azide_to_amine_reaction = checker.check_reaction(
                "Azide to amine reduction (Staudinger)", rsmi
            )

            # Detect amine to azide conversion
            if (reactants_have_amine and product_has_azide) or amine_to_azide_reaction:
                amine_to_azide_conversion = True
                print(f"Detected amine to azide conversion in reaction: {rsmi}")

            # Detect azide to amine conversion
            if (reactants_have_azide and product_has_amine) or azide_to_amine_reaction:
                azide_to_amine_conversion = True
                print(f"Detected azide to amine conversion in reaction: {rsmi}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is detected if either conversion is found
    # Prioritize amine-to-azide as it's the primary focus of the function
    result = amine_to_azide_conversion or azide_to_amine_conversion
    print(f"Amine-azide interconversion strategy detected: {result}")
    return result
