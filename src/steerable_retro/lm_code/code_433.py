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
    This function detects a sequential functional group transformation strategy:
    aldehyde → oxime → nitrile in a linear synthesis.
    """
    # Track the sequence of transformations
    transformation_sequence = []

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        current_path = path + [node]

        # Process molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            print(f"Depth {depth}, Processing molecule: {mol_smiles}")

            # Check for functional groups
            is_aldehyde = checker.check_fg("Aldehyde", mol_smiles)
            is_oxime = checker.check_fg("Oxime", mol_smiles)
            is_nitrile = checker.check_fg("Nitrile", mol_smiles)

            if is_aldehyde:
                print(f"Found aldehyde at depth {depth}: {mol_smiles}")
            if is_oxime:
                print(f"Found oxime at depth {depth}: {mol_smiles}")
            if is_nitrile:
                print(f"Found nitrile at depth {depth}: {mol_smiles}")

            # Store molecule with its functional groups for later analysis
            if is_aldehyde or is_oxime or is_nitrile:
                transformation_sequence.append(
                    {
                        "depth": depth,
                        "smiles": mol_smiles,
                        "is_aldehyde": is_aldehyde,
                        "is_oxime": is_oxime,
                        "is_nitrile": is_nitrile,
                        "node": node,
                    }
                )

        # Process reaction nodes to check reaction types
        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rxn_smiles = node["metadata"]["rsmi"]
            print(f"Depth {depth}, Processing reaction: {rxn_smiles}")

            # Check if this is an aldehyde to oxime reaction
            aldehyde_to_oxime = checker.check_reaction("Ketone/aldehyde to hydrazone", rxn_smiles)
            if aldehyde_to_oxime:
                print(f"Found aldehyde to oxime reaction at depth {depth}")

            # Check if this is an oxime to nitrile reaction
            # This could be a dehydration reaction
            oxime_to_nitrile = False

            # Try to identify reactants and products
            try:
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                # Check if reactants contain oxime and product contains nitrile
                reactant_has_oxime = any(checker.check_fg("Oxime", r) for r in reactants)
                product_has_nitrile = checker.check_fg("Nitrile", product)

                if reactant_has_oxime and product_has_nitrile:
                    oxime_to_nitrile = True
                    print(f"Found oxime to nitrile reaction at depth {depth}")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth to analyze the sequence (lower depth = later in synthesis)
    transformation_sequence.sort(key=lambda x: x["depth"])

    # Check if we have the sequence: aldehyde → oxime → nitrile
    # The sequence should be in reverse order since we're traversing retrosynthetically
    has_sequence = False

    if len(transformation_sequence) >= 3:
        print("Analyzing transformation sequence:")
        for item in transformation_sequence:
            print(
                f"Depth {item['depth']}: Aldehyde={item['is_aldehyde']}, Oxime={item['is_oxime']}, Nitrile={item['is_nitrile']}, SMILES={item['smiles']}"
            )

        # Look for the pattern: nitrile (late stage) → oxime (middle) → aldehyde (early stage)
        for i in range(len(transformation_sequence) - 2):
            if (
                transformation_sequence[i]["is_nitrile"]
                and transformation_sequence[i + 1]["is_oxime"]
                and transformation_sequence[i + 2]["is_aldehyde"]
            ):
                has_sequence = True
                print(
                    f"Found sequence: Nitrile (depth {transformation_sequence[i]['depth']}) → "
                    f"Oxime (depth {transformation_sequence[i+1]['depth']}) → "
                    f"Aldehyde (depth {transformation_sequence[i+2]['depth']})"
                )
                break

    print(f"Sequential aldehyde → oxime → nitrile strategy detected: {has_sequence}")
    return has_sequence
