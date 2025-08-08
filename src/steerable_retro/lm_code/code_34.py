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
    This function detects transformation of an aldehyde to a terminal alkyne.
    Typically a late-stage functional group interconversion.
    """
    transformation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal transformation_detected

        # Consider reactions at depth < 3 as late-stage
        is_late_stage = depth < 3

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde in reactants using checker function
                has_aldehyde = any(checker.check_fg("Aldehyde", reactant) for reactant in reactants)

                # Check for terminal alkyne in product using checker function
                has_terminal_alkyne = checker.check_fg("Alkyne", product)

                # Check if aldehyde is gone in product
                product_has_aldehyde = checker.check_fg("Aldehyde", product)

                # Check if alkyne already exists in reactants
                reactants_have_alkyne = any(
                    checker.check_fg("Alkyne", reactant) for reactant in reactants
                )

                # Check if this is a relevant reaction type
                is_homologation = checker.check_reaction(
                    "Homologation of aldehydes with formaldehyde", rsmi
                )

                # Verify transformation: aldehyde should be consumed, alkyne should be new
                if (
                    has_aldehyde
                    and has_terminal_alkyne
                    and is_late_stage
                    and not reactants_have_alkyne
                ):
                    if is_homologation:
                        print(
                            f"Detected aldehyde to terminal alkyne transformation via homologation at depth {depth}"
                        )
                        transformation_detected = True
                    elif not product_has_aldehyde:
                        # If specific reaction check fails, look for general pattern
                        # The aldehyde must be consumed and the alkyne must be new
                        print(
                            f"Detected general aldehyde to terminal alkyne transformation at depth {depth}"
                        )
                        transformation_detected = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return transformation_detected
