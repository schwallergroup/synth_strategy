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
    This function detects the reduction of a ketone to an amine in the synthesis.
    """
    has_ketone_reduction = False

    def dfs_traverse(node):
        nonlocal has_ketone_reduction

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ketone in reactants
                ketone_in_reactants = False
                for reactant in reactants:
                    if checker.check_fg("Ketone", reactant):
                        ketone_in_reactants = True
                        break

                # Check for amine in product
                amine_in_product = checker.check_fg("Primary amine", product) or checker.check_fg(
                    "Secondary amine", product
                )

                # Check for appropriate reaction types
                is_reductive_amination = checker.check_reaction(
                    "Reductive amination with ketone", rsmi
                )
                is_reduction_reaction = (
                    checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                    and amine_in_product
                )

                # Additional check for other reduction pathways
                is_other_reduction = False
                if ketone_in_reactants and amine_in_product:
                    # Check if this is a reduction reaction by looking for common patterns
                    if "NaBH" in rsmi or "LiAlH" in rsmi or "H2" in rsmi or "NH3" in rsmi:
                        is_other_reduction = True

                if is_reductive_amination or is_reduction_reaction or is_other_reduction:
                    has_ketone_reduction = True
                    print(f"Detected ketone to amine reduction: {rsmi}")
                    print(f"Reductive amination: {is_reductive_amination}")
                    print(f"Reduction reaction: {is_reduction_reaction}")
                    print(f"Other reduction: {is_other_reduction}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_ketone_reduction
