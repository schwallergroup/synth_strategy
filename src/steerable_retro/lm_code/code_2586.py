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
    Detects linear synthesis of benzylic nitrile through sequential functional group
    transformations: aldehyde → alcohol → chloride → nitrile
    """
    # Track transformations in sequence
    transformations = []

    def is_benzylic(smiles, fg_type=None):
        """Check if a functional group is in benzylic position"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Different SMARTS patterns based on functional group type
        if fg_type == "Aldehyde":
            return mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH]=O"))
        elif fg_type == "Primary alcohol":
            return mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH2]O"))
        elif fg_type == "Primary halide":
            return (
                mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH2]Cl"))
                or mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH2]Br"))
                or mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH2]I"))
            )
        elif fg_type == "Nitrile":
            return mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH2]C#N"))
        else:
            # Generic benzylic check
            return mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH2]"))

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Depth {depth} - Examining reaction: {rsmi}")

                # Check for nitrile formation from chloride
                if (
                    any(
                        checker.check_fg("Primary halide", r) and is_benzylic(r, "Primary halide")
                        for r in reactants
                    )
                    and checker.check_fg("Nitrile", product)
                    and is_benzylic(product, "Nitrile")
                ):
                    transformations.append("chloride_to_nitrile")
                    print(f"Depth {depth} - Detected: chloride → nitrile")

                # Check for chloride formation from alcohol
                if (
                    any(
                        checker.check_fg("Primary alcohol", r) and is_benzylic(r, "Primary alcohol")
                        for r in reactants
                    )
                    and checker.check_fg("Primary halide", product)
                    and is_benzylic(product, "Primary halide")
                ):
                    transformations.append("alcohol_to_chloride")
                    print(f"Depth {depth} - Detected: alcohol → chloride")

                # Check for alcohol formation from aldehyde
                if (
                    any(
                        checker.check_fg("Aldehyde", r) and is_benzylic(r, "Aldehyde")
                        for r in reactants
                    )
                    and checker.check_fg("Primary alcohol", product)
                    and is_benzylic(product, "Primary alcohol")
                ):
                    transformations.append("aldehyde_to_alcohol")
                    print(f"Depth {depth} - Detected: aldehyde → alcohol")

                # Additional check for aldehyde to alcohol transformation
                # Sometimes the reaction might be detected in the opposite direction due to retrosynthetic analysis
                if (
                    checker.check_fg("Aldehyde", product)
                    and is_benzylic(product, "Aldehyde")
                    and any(
                        checker.check_fg("Primary alcohol", r) and is_benzylic(r, "Primary alcohol")
                        for r in reactants
                    )
                ):
                    if checker.check_reaction(
                        "Oxidation of aldehydes to carboxylic acids", rsmi
                    ) or checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                    ):
                        transformations.append("aldehyde_to_alcohol")
                        print(f"Depth {depth} - Detected (reverse): aldehyde → alcohol")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Detected transformations: {transformations}")

    # Check if all three transformations are present
    required_transformations = ["chloride_to_nitrile", "alcohol_to_chloride", "aldehyde_to_alcohol"]
    sequence_present = all(t in transformations for t in required_transformations)

    print(f"Linear benzylic nitrile synthesis detected: {sequence_present}")

    return sequence_present
