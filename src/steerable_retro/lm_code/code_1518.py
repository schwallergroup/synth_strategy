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
    This function detects the use of ester hydrolysis to generate carboxylic acids.
    """
    # Track if we found an ester hydrolysis
    found_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal found_ester_hydrolysis

        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction: {rsmi}")

                # Extract reactants and product
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # First check predefined reaction types
                is_hydrolysis_reaction = (
                    checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                )

                if is_hydrolysis_reaction:
                    print(f"Reaction matches a known hydrolysis type")

                # Check for ester/thioester in reactants
                has_ester_reactant = any(checker.check_fg("Ester", r) for r in reactants_smiles)
                has_thioester_reactant = any(
                    checker.check_fg("Carbo-thioester", r) for r in reactants_smiles
                )

                # Check for carboxylic acid in product
                has_carboxylic_acid_product = checker.check_fg("Carboxylic acid", product_smiles)

                # Count esters in reactants and product
                ester_count_reactants = sum(
                    1 for r in reactants_smiles if checker.check_fg("Ester", r)
                )
                ester_count_product = 1 if checker.check_fg("Ester", product_smiles) else 0

                print(
                    f"Reactants - Ester: {has_ester_reactant}, Thioester: {has_thioester_reactant}, Ester count: {ester_count_reactants}"
                )
                print(
                    f"Product - Carboxylic acid: {has_carboxylic_acid_product}, Ester count: {ester_count_product}"
                )

                # Check if this is an ester hydrolysis:
                # 1. Either a known hydrolysis reaction type, OR
                # 2. Has ester/thioester in reactants, carboxylic acid in product, and fewer esters in product
                if is_hydrolysis_reaction or (
                    (has_ester_reactant or has_thioester_reactant)
                    and has_carboxylic_acid_product
                    and ester_count_product < ester_count_reactants
                ):
                    print(f"Found ester hydrolysis reaction: {rsmi}")
                    found_ester_hydrolysis = True
                else:
                    print(f"Not an ester hydrolysis reaction")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    print("Starting traversal of synthetic route")
    dfs_traverse(route)
    print(f"Ester hydrolysis found: {found_ester_hydrolysis}")

    return found_ester_hydrolysis
