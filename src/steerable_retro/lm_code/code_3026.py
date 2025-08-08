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
    This function detects a strategy for synthesizing phosphonate-containing aniline compounds.
    """
    # Track if we found the target structure and the synthesis strategy
    has_target_structure = False
    has_phosphonate_strategy = False
    has_aniline_strategy = False

    # Track the depth of each node to identify the final product
    node_depths = {}

    # First pass: calculate depths
    def calculate_depths(node, depth=0):
        node_id = id(node)
        node_depths[node_id] = depth

        for child in node.get("children", []):
            calculate_depths(child, depth + 1)

    calculate_depths(route)

    # Find the root node (final product) - it has the minimum depth
    min_depth = min(node_depths.values()) if node_depths else 0

    # Helper function to check for phosphonate groups
    def has_phosphonate_group(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Check for any phosphorus atom
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "P":
                        return True
                # Also check for phosphate ester specifically
                return checker.check_fg("Phosphate ester", smiles)
            return False
        except:
            return False

    # Second pass: analyze the synthesis route
    def analyze_route(node):
        nonlocal has_target_structure, has_phosphonate_strategy, has_aniline_strategy

        if node["type"] == "mol" and "smiles" in node:
            node_id = id(node)
            mol_smiles = node["smiles"]

            # Check if this is the final product (root node or close to it)
            if node_depths[node_id] <= min_depth + 1:
                # Check for phosphonate group
                has_phosphonate = has_phosphonate_group(mol_smiles)
                # Check for aniline group
                has_aniline = checker.check_fg("Aniline", mol_smiles)

                if has_phosphonate and has_aniline:
                    has_target_structure = True
                    print(
                        f"Found target structure with phosphonate and aniline groups: {mol_smiles}"
                    )

        # Check reactions for synthesis strategy
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction introduces or modifies phosphonate group
            product_has_phosphonate = has_phosphonate_group(product)
            reactants_have_phosphonate = any(has_phosphonate_group(r) for r in reactants)

            # Check if this reaction introduces or modifies aniline group
            product_has_aniline = checker.check_fg("Aniline", product)
            reactants_have_aniline = any(checker.check_fg("Aniline", r) for r in reactants)
            reactants_have_nitro = any(checker.check_fg("Nitro group", r) for r in reactants)

            # Check for phosphonate formation or modification
            if product_has_phosphonate:
                if not reactants_have_phosphonate:
                    print(f"Found phosphonate formation reaction: {rsmi}")
                    has_phosphonate_strategy = True
                elif reactants_have_phosphonate:
                    # Check if this is a modification of an existing phosphonate
                    print(f"Found phosphonate modification reaction: {rsmi}")
                    has_phosphonate_strategy = True

            # Check for aniline formation or modification
            if product_has_aniline:
                if not reactants_have_aniline:
                    print(f"Found aniline formation reaction: {rsmi}")
                    has_aniline_strategy = True

                    # Specifically check for nitro reduction to aniline
                    if reactants_have_nitro and checker.check_reaction(
                        "Reduction of nitro groups to amines", rsmi
                    ):
                        print(f"Found nitro reduction to aniline: {rsmi}")
                        has_aniline_strategy = True
                    # Check for other amine formation reactions
                    elif (
                        checker.check_reaction("Reductive amination with aldehyde", rsmi)
                        or checker.check_reaction("Reductive amination with ketone", rsmi)
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )
                    ):
                        print(f"Found aniline formation through alternative reaction: {rsmi}")
                        has_aniline_strategy = True
                elif reactants_have_aniline:
                    # Check if this is a modification of an existing aniline
                    print(f"Found aniline modification reaction: {rsmi}")
                    has_aniline_strategy = True

        # Traverse children
        for child in node.get("children", []):
            analyze_route(child)

    # Start traversal from the root
    analyze_route(route)

    # Return True if both target structure and synthesis strategies are found
    has_synthesis_strategy = has_phosphonate_strategy or has_aniline_strategy
    result = has_target_structure or has_synthesis_strategy
    print(f"Target structure found: {has_target_structure}")
    print(f"Phosphonate strategy found: {has_phosphonate_strategy}")
    print(f"Aniline strategy found: {has_aniline_strategy}")
    print(f"Final result: {result}")

    return result
