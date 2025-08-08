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
    This function detects a synthetic strategy involving construction of a tertiary carbon
    center via nucleophilic addition to a carbonyl compound.
    """
    # Initialize tracking variables
    has_tertiary_carbon_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_tertiary_carbon_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for Grignard or organolithium reactions
            if (
                checker.check_reaction("Grignard from aldehyde to alcohol", rsmi)
                or checker.check_reaction("Grignard from ketone to alcohol", rsmi)
                or checker.check_reaction(
                    "Reaction of alkyl halides with organometallic coumpounds", rsmi
                )
            ):

                # Verify tertiary alcohol formation in product
                if checker.check_fg("Tertiary alcohol", product_part):
                    # Check that the tertiary alcohol wasn't already present in reactants
                    if not checker.check_fg("Tertiary alcohol", reactants_part):
                        print(
                            f"Detected tertiary carbon formation via carbonyl addition at depth {depth}"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        has_tertiary_carbon_formation = True

            # Check for other carbonyl addition reactions that might form tertiary carbons
            elif (
                checker.check_fg("Aldehyde", reactants_part)
                or checker.check_fg("Ketone", reactants_part)
                or checker.check_fg("Ester", reactants_part)
                or checker.check_fg("Carboxylic acid", reactants_part)
            ):

                # Verify tertiary alcohol formation in product
                if checker.check_fg("Tertiary alcohol", product_part):
                    # Check that the tertiary alcohol wasn't already present in reactants
                    if not checker.check_fg("Tertiary alcohol", reactants_part):
                        print(
                            f"Detected tertiary carbon formation via carbonyl addition at depth {depth}"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        has_tertiary_carbon_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_tertiary_carbon_formation
