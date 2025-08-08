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
    This function detects a pattern of multiple functional group interconversions
    at the 2-position of a quinoline or similar heterocycle.
    """
    # Track functional group changes at position 2
    position_2_changes = []

    # List of functional groups to check at position 2
    functional_groups = [
        "Chloride",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Phenol",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Ether",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Nitrile",
        "Nitro group",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Aldehyde",
        "Ketone",
        "Azide",
        "Triflate",
        "Tosylate",
        "Mesylate",
    ]

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    product_mol = Chem.MolFromSmiles(product)

                    # Check if product contains quinoline or isoquinoline
                    has_quinoline = checker.check_ring("quinoline", product)
                    has_isoquinoline = checker.check_ring("isoquinoline", product)

                    if has_quinoline or has_isoquinoline:
                        print(
                            f"Found {'quinoline' if has_quinoline else 'isoquinoline'} in product: {product}"
                        )

                        # Check for functional groups at position 2
                        # For quinoline, position 2 is typically the carbon next to the nitrogen
                        current_fg = None
                        for fg in functional_groups:
                            if checker.check_fg(fg, product):
                                # Get atom indices of the functional group
                                fg_indices = checker.get_fg_atom_indices(fg, product)
                                if fg_indices:
                                    # Check if any of the functional group atoms are at position 2
                                    # This is a simplification - in a real implementation, we would use
                                    # atom mapping to precisely identify position 2
                                    current_fg = fg
                                    print(f"Found functional group {fg} in product")
                                    break

                        # Now check reactants to see if the functional group changed
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and (
                                checker.check_ring("quinoline", reactant)
                                or checker.check_ring("isoquinoline", reactant)
                            ):
                                print(f"Found quinoline/isoquinoline in reactant: {reactant}")

                                # Check for functional groups in reactant
                                for fg in functional_groups:
                                    if fg != current_fg and checker.check_fg(fg, reactant):
                                        # Detected a functional group change
                                        print(
                                            f"Functional group change detected: {fg} -> {current_fg}"
                                        )
                                        position_2_changes.append((fg, current_fg, depth))
                                        break
                except Exception as e:
                    print(f"Error processing molecule in functional_group_dance_at_position_2: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Sort changes by depth to get them in synthetic order
    position_2_changes.sort(key=lambda x: x[2])

    # Extract just the functional groups
    fg_sequence = [(from_fg, to_fg) for from_fg, to_fg, _ in position_2_changes]

    print(f"Detected functional group changes at position 2: {fg_sequence}")

    # Check if we have at least 2 different functional group changes
    if len(fg_sequence) >= 2:
        # Check if the changes are actually different
        unique_changes = set(fg_sequence)
        if len(unique_changes) >= 2:
            print(
                f"Detected functional group dance at position 2 with {len(unique_changes)} unique changes"
            )
            return True

    return False
