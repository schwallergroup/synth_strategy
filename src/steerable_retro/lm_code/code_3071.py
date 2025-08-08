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
    This function detects a strategy involving ester hydrolysis.
    """
    has_ester_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal has_ester_hydrolysis

        if node["type"] == "reaction":
            # Extract product and reactants
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is an ester hydrolysis reaction - check multiple reaction types
                is_hydrolysis = (
                    checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                )

                if is_hydrolysis:
                    print(f"Found potential ester hydrolysis reaction at depth {depth}")

                    # Verify ester and carboxylic acid presence
                    has_ester = False
                    has_acid = False

                    # Check reactants
                    for reactant in reactants_smiles:
                        if checker.check_fg("Ester", reactant):
                            has_ester = True
                            print(f"Found ester in reactant: {reactant}")
                        if checker.check_fg("Carboxylic acid", reactant):
                            has_acid = True
                            print(f"Found carboxylic acid in reactant: {reactant}")

                    # Check product
                    if checker.check_fg("Ester", product_smiles):
                        has_ester = True
                        print(f"Found ester in product: {product_smiles}")
                    if checker.check_fg("Carboxylic acid", product_smiles):
                        has_acid = True
                        print(f"Found carboxylic acid in product: {product_smiles}")

                    # Confirm ester hydrolysis if both functional groups are present
                    if has_ester and has_acid:
                        has_ester_hydrolysis = True
                        print(f"Confirmed ester hydrolysis at depth {depth}")

                # Even if not detected by reaction type, check for pattern of ester to acid conversion
                else:
                    # Check if reactants contain ester and product contains acid
                    has_ester_in_reactants = any(
                        checker.check_fg("Ester", r) for r in reactants_smiles
                    )
                    has_acid_in_product = checker.check_fg("Carboxylic acid", product_smiles)

                    # Check if product contains ester and reactants contain acid (retrosynthetic direction)
                    has_ester_in_product = checker.check_fg("Ester", product_smiles)
                    has_acid_in_reactants = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )

                    # Check for ester hydrolysis pattern in either direction
                    if (has_ester_in_reactants and has_acid_in_product) or (
                        has_ester_in_product and has_acid_in_reactants
                    ):
                        print(
                            f"Detected ester hydrolysis pattern at depth {depth} without specific reaction type"
                        )
                        has_ester_hydrolysis = True

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: has_ester_hydrolysis = {has_ester_hydrolysis}")

    return has_ester_hydrolysis
