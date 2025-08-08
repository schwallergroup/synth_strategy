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
    This function detects a synthetic strategy involving sequential functional group
    interconversions on the same carbon atom, specifically acid chloride → alcohol →
    aldehyde → amine, while preserving a heterocyclic scaffold.
    """
    # Track the sequence of transformations
    transformations = []

    # Track if we have a heterocyclic scaffold
    has_heterocyclic_scaffold = False

    # Define heterocyclic rings to check
    heterocyclic_rings = [
        "benzothiophene",
        "indole",
        "benzoxazole",
        "benzimidazole",
        "quinoline",
        "isoquinoline",
        "benzofuran",
        "carbazole",
        "purine",
        "thiazole",
        "oxazole",
        "imidazole",
        "pyridine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocyclic_scaffold

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for heterocyclic scaffold
            if not has_heterocyclic_scaffold:
                for ring in heterocyclic_rings:
                    if checker.check_ring(ring, mol_smiles):
                        print(f"Found heterocyclic scaffold: {ring}")
                        has_heterocyclic_scaffold = True
                        break

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acid chloride to alcohol transformation
            if any(
                checker.check_fg("Acyl halide", reactant) for reactant in reactants
            ) and checker.check_fg("Primary alcohol", product):
                print(f"Found transformation: Acyl halide → Primary alcohol")
                transformations.append(("acid_chloride_to_alcohol", depth))

            # Check for alcohol to aldehyde transformation
            if any(
                checker.check_fg("Primary alcohol", reactant) for reactant in reactants
            ) and checker.check_fg("Aldehyde", product):
                print(f"Found transformation: Primary alcohol → Aldehyde")
                transformations.append(("alcohol_to_aldehyde", depth))

            # Check for aldehyde to amine transformation (reductive amination)
            if checker.check_reaction("Reductive amination with aldehyde", rsmi) or (
                any(checker.check_fg("Aldehyde", reactant) for reactant in reactants)
                and any(checker.check_fg("Primary amine", reactant) for reactant in reactants)
                and (
                    checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )
            ):
                print(f"Found transformation: Aldehyde → Amine (Reductive amination)")
                transformations.append(("aldehyde_to_amine", depth))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth (to get the correct sequence)
    transformations.sort(key=lambda x: x[1], reverse=True)

    # Extract just the transformation types in sequence
    transformation_sequence = [t[0] for t in transformations]
    print(f"Transformation sequence: {transformation_sequence}")

    # Check if we have the correct sequence
    has_correct_sequence = False
    if len(transformation_sequence) >= 3:
        for i in range(len(transformation_sequence) - 2):
            if (
                transformation_sequence[i] == "acid_chloride_to_alcohol"
                and transformation_sequence[i + 1] == "alcohol_to_aldehyde"
                and transformation_sequence[i + 2] == "aldehyde_to_amine"
            ):
                has_correct_sequence = True
                break

    # Check if we found the complete strategy
    strategy_found = has_correct_sequence and has_heterocyclic_scaffold

    print(f"Sequential functional group interconversion strategy detected: {strategy_found}")
    print(f"  - Correct transformation sequence: {has_correct_sequence}")
    print(f"  - Heterocyclic scaffold found: {has_heterocyclic_scaffold}")

    return strategy_found
