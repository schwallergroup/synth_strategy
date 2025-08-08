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
    This function detects a synthetic strategy involving an ester to aldehyde
    conversion as part of a functional group transformation sequence.
    """
    has_ester_to_aldehyde = False

    # Track molecules with esters and aldehydes at each depth
    molecules_with_ester = {}
    molecules_with_aldehyde = {}

    def dfs_traverse(node, depth=0, path=None):
        nonlocal has_ester_to_aldehyde

        if path is None:
            path = []

        # For molecule nodes, check and track functional groups
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Track molecules with functional groups
            if checker.check_fg("Ester", mol_smiles):
                if depth not in molecules_with_ester:
                    molecules_with_ester[depth] = []
                molecules_with_ester[depth].append(mol_smiles)

            if checker.check_fg("Aldehyde", mol_smiles):
                if depth not in molecules_with_aldehyde:
                    molecules_with_aldehyde[depth] = []
                molecules_with_aldehyde[depth].append(mol_smiles)

        # For reaction nodes, check for ester to aldehyde transformations
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Split reactants for individual checking
            reactant_list = reactants_smiles.split(".")

            # Check for functional groups in reactants and products
            has_ester_in_reactants = any(checker.check_fg("Ester", r) for r in reactant_list)
            has_aldehyde_in_reactants = any(checker.check_fg("Aldehyde", r) for r in reactant_list)
            has_ester_in_product = checker.check_fg("Ester", product_smiles)
            has_aldehyde_in_product = checker.check_fg("Aldehyde", product_smiles)

            # Check for relevant reactions that could be part of ester to aldehyde conversion
            # Direct ester to aldehyde conversion
            if has_ester_in_reactants and has_aldehyde_in_product:
                print(f"Detected direct ester to aldehyde conversion at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                has_ester_to_aldehyde = True

            # Check for specific reactions that could be part of the pathway
            # Reduction of ester to aldehyde
            if checker.check_reaction("Reduction of ester to primary alcohol", rsmi):
                print(f"Detected ester reduction to alcohol at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                # This could be part of a multi-step pathway
                if depth + 1 in molecules_with_aldehyde:
                    print(f"Found potential multi-step ester to aldehyde pathway")
                    has_ester_to_aldehyde = True

            # Oxidation of alcohol to aldehyde
            if (
                checker.check_reaction(
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                )
                and has_aldehyde_in_product
            ):
                print(f"Detected alcohol oxidation to aldehyde at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                # Check if there's an ester at a higher depth (earlier in synthesis)
                for d in molecules_with_ester:
                    if d > depth:
                        print(f"Found potential multi-step ester to aldehyde pathway")
                        has_ester_to_aldehyde = True
                        break

            # Reduction of carboxylic acid to aldehyde (via alcohol)
            if checker.check_reaction("Reduction of carboxylic acid to primary alcohol", rsmi):
                print(f"Detected carboxylic acid reduction to alcohol at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                # This could be part of a multi-step pathway
                if depth + 1 in molecules_with_aldehyde:
                    print(f"Found potential multi-step acid to aldehyde pathway")
                    has_ester_to_aldehyde = True

            # Ester hydrolysis to carboxylic acid
            if (
                checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                )
                and has_ester_in_reactants
            ):
                print(f"Detected ester hydrolysis to carboxylic acid at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                # Check if there's an aldehyde at a lower depth (later in synthesis)
                for d in molecules_with_aldehyde:
                    if d < depth:
                        print(f"Found potential multi-step ester to aldehyde pathway")
                        has_ester_to_aldehyde = True
                        break

            # Aldehyde oxidation to ester (in retrosynthesis: ester reduction to aldehyde)
            if (
                checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                and has_aldehyde_in_reactants
            ):
                print(f"Detected aldehyde oxidation to carboxylic acid at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                has_ester_to_aldehyde = True

            # Esterification (in retrosynthesis: ester to carboxylic acid)
            if (
                checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                and has_ester_in_product
            ):
                print(f"Detected esterification at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                # Check if there's an aldehyde at a lower depth (later in synthesis)
                for d in molecules_with_aldehyde:
                    if d < depth:
                        print(f"Found potential multi-step ester to aldehyde pathway")
                        has_ester_to_aldehyde = True
                        break

        # Add current node to path and traverse children
        path.append(node)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, path)
        path.pop()

    # Start traversal from the root
    dfs_traverse(route)

    # Check for multi-step pathways across the entire route
    if not has_ester_to_aldehyde:
        # If we have both esters and aldehydes in the route, check if they could be connected
        if molecules_with_ester and molecules_with_aldehyde:
            # Find the highest depth (earliest in synthesis) with an ester
            max_ester_depth = max(molecules_with_ester.keys())
            # Find the lowest depth (latest in synthesis) with an aldehyde
            min_aldehyde_depth = min(molecules_with_aldehyde.keys())

            # If the ester appears before the aldehyde in the synthesis route,
            # it could be part of an ester to aldehyde conversion strategy
            if max_ester_depth > min_aldehyde_depth:
                print(f"Detected potential multi-step ester to aldehyde pathway")
                print(f"Esters found at depths: {list(molecules_with_ester.keys())}")
                print(f"Aldehydes found at depths: {list(molecules_with_aldehyde.keys())}")
                has_ester_to_aldehyde = True

    return has_ester_to_aldehyde
