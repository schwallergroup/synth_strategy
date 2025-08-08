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
    This function detects a strategy involving thiocarbamate formation from an amine.
    In retrosynthesis, this means finding reactions where a thiocarbamate is converted back to an amine.
    """
    thiocarbamate_formation_found = False

    def dfs_traverse(node):
        nonlocal thiocarbamate_formation_found

        if node["type"] == "reaction" and not thiocarbamate_formation_found:
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis: reactants have thiocarbamate, product has amine
                # Thiocarbamate is a type of thioamide with O-C(=S)-N structure
                reactant_has_thioamide = any(checker.check_fg("Thioamide", r) for r in reactants)

                # Check if product contains a primary amine
                product_has_amine = checker.check_fg("Primary amine", product) or checker.check_fg(
                    "Secondary amine", product
                )

                # Additional check for thiourea reaction which can form thiocarbamates
                is_thiourea_reaction = checker.check_reaction("Thiourea", rsmi)

                if reactant_has_thioamide and product_has_amine:
                    print(f"Found thiocarbamate formation from amine in reaction: {rsmi}")
                    thiocarbamate_formation_found = True
                elif is_thiourea_reaction:
                    print(f"Found thiourea reaction which can form thiocarbamates: {rsmi}")
                    thiocarbamate_formation_found = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal if strategy not found yet
        if not thiocarbamate_formation_found:
            for child in node.get("children", []):
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return thiocarbamate_formation_found
