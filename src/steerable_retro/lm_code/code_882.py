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
    This function detects if the synthetic route uses azide chemistry as a key intermediate.
    It looks for azide functional groups and azide formation reactions.
    """
    azide_found = False

    # List of azide-related reactions to check
    azide_formation_reactions = [
        "Formation of Azides from halogens",
        "Formation of Azides from boronic acids",
        "Alcohol to azide",
        "Amine to azide",
    ]

    azide_utilizing_reactions = [
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Azide to amine reduction (Staudinger)",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
    ]

    def dfs_traverse(node):
        nonlocal azide_found

        if azide_found:
            return  # Early return if we already found azide chemistry

        if node["type"] == "mol":
            # Check if molecule contains azide group
            if checker.check_fg("Azide", node["smiles"]):
                azide_found = True
                print(f"Found azide functional group in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for azide formation reactions
                for reaction_name in azide_formation_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        azide_found = True
                        print(f"Found azide formation reaction ({reaction_name}): {rsmi}")
                        break

                # Check for reactions that utilize azides
                if not azide_found:
                    for reaction_name in azide_utilizing_reactions:
                        if checker.check_reaction(reaction_name, rsmi):
                            azide_found = True
                            print(f"Found azide-utilizing reaction ({reaction_name}): {rsmi}")
                            break

                # Additional check: if reactants or products contain azide groups
                if not azide_found:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant contains azide
                        reactant_has_azide = any(
                            checker.check_fg("Azide", reactant) for reactant in reactants
                        )

                        # Check if product contains azide
                        product_has_azide = checker.check_fg("Azide", product)

                        # If azide appears or disappears, it's azide chemistry
                        if reactant_has_azide or product_has_azide:
                            azide_found = True
                            print(f"Found reaction involving azide groups: {rsmi}")
                    except Exception as e:
                        print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return azide_found
