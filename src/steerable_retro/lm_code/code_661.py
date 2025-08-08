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
    This function detects if the synthesis involves the formation of an acid chloride
    as an intermediate step.
    """
    found_acid_chloride_formation = False

    def dfs_traverse(node):
        nonlocal found_acid_chloride_formation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]

            # Extract reactants and product
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]
            except:
                print(f"Error extracting reactants and product from: {rsmi}")
                return

            # Check if this is an acid chloride formation reaction using checker functions
            if (
                checker.check_reaction("Acyl chlorides from alcohols", rsmi)
                or checker.check_reaction("Alcohol to chloride_SOCl2", rsmi)
                or checker.check_reaction("Alcohol to chloride_PCl5_ortho", rsmi)
                or checker.check_reaction("Alcohol to chloride_POCl3", rsmi)
            ):

                # Verify that an acid chloride is formed by checking functional groups
                if checker.check_fg("Acyl halide", product) and any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants
                ):
                    found_acid_chloride_formation = True
                    print(f"Found acid chloride formation in reaction: {rsmi}")

            # Check for other reactions that might form acid chlorides
            elif any(
                checker.check_fg("Carboxylic acid", r) for r in reactants
            ) and checker.check_fg("Acyl halide", product):
                found_acid_chloride_formation = True
                print(f"Found acid chloride formation (general case) in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Acid chloride formation strategy: {found_acid_chloride_formation}")
    return found_acid_chloride_formation
