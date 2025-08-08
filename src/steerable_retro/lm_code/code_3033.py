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
    This function detects reductive amination pattern in the synthesis.
    """
    found_reductive_amination = False

    def dfs_traverse(node):
        nonlocal found_reductive_amination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a reductive amination reaction
            if (
                checker.check_reaction("Reductive amination with aldehyde", rsmi)
                or checker.check_reaction("Reductive amination with ketone", rsmi)
                or checker.check_reaction("Reductive amination with alcohol", rsmi)
            ):
                print(f"Found reductive amination reaction: {rsmi}")
                found_reductive_amination = True
            else:
                # Alternative check if reaction type check fails
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde and amine components
                has_aldehyde = any(checker.check_fg("Aldehyde", reactant) for reactant in reactants)
                has_ketone = any(checker.check_fg("Ketone", reactant) for reactant in reactants)
                has_amine = any(
                    checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                    for reactant in reactants
                )

                # Check for new C-N bond in product
                if (has_aldehyde or has_ketone) and has_amine:
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Check if product has a new C-N bond characteristic of reductive amination
                            if checker.check_fg("Secondary amine", product) or checker.check_fg(
                                "Tertiary amine", product
                            ):
                                print(f"Found likely reductive amination pattern: {rsmi}")
                                found_reductive_amination = True
                    except Exception as e:
                        print(f"Error analyzing product: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_reductive_amination
