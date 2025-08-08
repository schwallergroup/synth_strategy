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
    This function detects a synthetic strategy where a carboxylic acid is synthesized
    from a primary alcohol via diazo and acyl chloride intermediates, while preserving
    an ester functional group.
    """
    # Initialize flags to track the presence of key intermediates and transformations
    has_primary_alcohol = False
    has_diazo = False
    has_acyl_chloride = False
    has_carboxylic_acid = False
    has_ester_throughout = True

    # Track reactions between key intermediates
    alcohol_to_diazo_rxn = False
    diazo_to_acyl_chloride_rxn = False
    acyl_chloride_to_carboxylic_acid_rxn = False

    # Track molecules and their depths
    primary_alcohol_mol = None
    diazo_mol = None
    acyl_chloride_mol = None
    carboxylic_acid_mol = None

    # Track the sequence of transformations
    transformation_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal has_primary_alcohol, has_diazo, has_acyl_chloride, has_carboxylic_acid
        nonlocal has_ester_throughout, alcohol_to_diazo_rxn, diazo_to_acyl_chloride_rxn
        nonlocal acyl_chloride_to_carboxylic_acid_rxn, primary_alcohol_mol, diazo_mol
        nonlocal acyl_chloride_mol, carboxylic_acid_mol

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for functional groups
            if checker.check_fg("Primary alcohol", mol_smiles) and checker.check_fg(
                "Ester", mol_smiles
            ):
                has_primary_alcohol = True
                primary_alcohol_mol = mol_smiles
                print(f"Found primary alcohol at depth {depth}: {mol_smiles}")

            if checker.check_fg("Diazo", mol_smiles) and checker.check_fg("Ester", mol_smiles):
                has_diazo = True
                diazo_mol = mol_smiles
                print(f"Found diazo with ester at depth {depth}: {mol_smiles}")

            if checker.check_fg("Acyl halide", mol_smiles) and checker.check_fg(
                "Ester", mol_smiles
            ):
                has_acyl_chloride = True
                acyl_chloride_mol = mol_smiles
                print(f"Found acyl chloride with ester at depth {depth}: {mol_smiles}")

            if checker.check_fg("Carboxylic acid", mol_smiles) and checker.check_fg(
                "Ester", mol_smiles
            ):
                has_carboxylic_acid = True
                carboxylic_acid_mol = mol_smiles
                print(f"Found carboxylic acid with ester at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]
            reactants = reactants_str.split(".")

            # Check for specific reactions
            # Alcohol to diazo reaction
            if checker.check_fg("Primary alcohol", product) and any(
                checker.check_fg("Diazo", r) for r in reactants
            ):
                if checker.check_reaction(
                    "O-alkylation of carboxylic acids with diazo compounds", rsmi
                ) or checker.check_reaction("O-alkylation of amides with diazo compounds", rsmi):
                    alcohol_to_diazo_rxn = True
                    transformation_sequence.append(("alcohol_to_diazo", depth))
                    print(f"Found alcohol to diazo reaction at depth {depth}: {rsmi}")

            # Diazo to acyl chloride reaction
            if checker.check_fg("Diazo", product) and any(
                checker.check_fg("Acyl halide", r) for r in reactants
            ):
                if (
                    checker.check_reaction("Acyl chlorides from alcohols", rsmi)
                    or checker.check_reaction("Acyl bromides from alcohols", rsmi)
                    or checker.check_reaction("Acyl iodides from alcohols", rsmi)
                ):
                    diazo_to_acyl_chloride_rxn = True
                    transformation_sequence.append(("diazo_to_acyl_chloride", depth))
                    print(f"Found diazo to acyl chloride reaction at depth {depth}: {rsmi}")

            # Acyl chloride to carboxylic acid reaction
            if checker.check_fg("Acyl halide", product) and any(
                checker.check_fg("Carboxylic acid", r) for r in reactants
            ):
                if (
                    checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                ):
                    acyl_chloride_to_carboxylic_acid_rxn = True
                    transformation_sequence.append(("acyl_chloride_to_carboxylic_acid", depth))
                    print(
                        f"Found acyl chloride to carboxylic acid reaction at depth {depth}: {rsmi}"
                    )

            # Manual check for reactions if reaction checkers don't catch them
            if (
                not alcohol_to_diazo_rxn
                and checker.check_fg("Primary alcohol", product)
                and any(checker.check_fg("Diazo", r) for r in reactants)
            ):
                alcohol_to_diazo_rxn = True
                transformation_sequence.append(("alcohol_to_diazo", depth))
                print(f"Found alcohol to diazo reaction (manual check) at depth {depth}: {rsmi}")

            if (
                not diazo_to_acyl_chloride_rxn
                and checker.check_fg("Diazo", product)
                and any(checker.check_fg("Acyl halide", r) for r in reactants)
            ):
                diazo_to_acyl_chloride_rxn = True
                transformation_sequence.append(("diazo_to_acyl_chloride", depth))
                print(
                    f"Found diazo to acyl chloride reaction (manual check) at depth {depth}: {rsmi}"
                )

            if (
                not acyl_chloride_to_carboxylic_acid_rxn
                and checker.check_fg("Acyl halide", product)
                and any(checker.check_fg("Carboxylic acid", r) for r in reactants)
            ):
                acyl_chloride_to_carboxylic_acid_rxn = True
                transformation_sequence.append(("acyl_chloride_to_carboxylic_acid", depth))
                print(
                    f"Found acyl chloride to carboxylic acid reaction (manual check) at depth {depth}: {rsmi}"
                )

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if all main pathway molecules have an ester group
    if primary_alcohol_mol and diazo_mol and acyl_chloride_mol and carboxylic_acid_mol:
        has_ester_throughout = (
            checker.check_fg("Ester", primary_alcohol_mol)
            and checker.check_fg("Ester", diazo_mol)
            and checker.check_fg("Ester", acyl_chloride_mol)
            and checker.check_fg("Ester", carboxylic_acid_mol)
        )

    # Check if the synthesis follows the expected pattern
    intermediates_present = (
        has_primary_alcohol and has_diazo and has_acyl_chloride and has_carboxylic_acid
    )
    reactions_present = (
        alcohol_to_diazo_rxn and diazo_to_acyl_chloride_rxn and acyl_chloride_to_carboxylic_acid_rxn
    )

    # Check if the transformation sequence is in the correct order
    correct_order = True
    if transformation_sequence:
        # Sort by depth (ascending)
        sorted_transformations = sorted(transformation_sequence, key=lambda x: x[1])
        # In retrosynthesis, we traverse from product to reactants
        expected_order = [
            "alcohol_to_diazo",
            "diazo_to_acyl_chloride",
            "acyl_chloride_to_carboxylic_acid",
        ]
        actual_order = [t[0] for t in sorted_transformations]

        print(f"Transformation sequence: {actual_order}")
        print(f"Expected order: {expected_order}")

        # Check if the actual order contains all expected transformations in the right order
        if len(actual_order) >= len(expected_order):
            for i, expected in enumerate(expected_order):
                if expected not in actual_order:
                    correct_order = False
                    print(f"Missing transformation: {expected}")
                    break
                if i > 0 and actual_order.index(expected) < actual_order.index(
                    expected_order[i - 1]
                ):
                    correct_order = False
                    print(f"Incorrect order: {expected} should come after {expected_order[i-1]}")
                    break
        else:
            correct_order = False
            print("Incomplete transformation sequence")
    else:
        correct_order = False
        print("No transformations detected")

    correct_sequence = (
        intermediates_present and reactions_present and has_ester_throughout and correct_order
    )

    if correct_sequence:
        print("Detected carboxylic acid synthesis via diazo and acyl chloride intermediates")
    else:
        print("Strategy not detected. Missing elements:")
        if not has_primary_alcohol:
            print("- No primary alcohol with ester")
        if not has_diazo:
            print("- No diazo intermediate with ester")
        if not has_acyl_chloride:
            print("- No acyl chloride intermediate with ester")
        if not has_carboxylic_acid:
            print("- No carboxylic acid product with ester")
        if not alcohol_to_diazo_rxn:
            print("- No alcohol to diazo reaction")
        if not diazo_to_acyl_chloride_rxn:
            print("- No diazo to acyl chloride reaction")
        if not acyl_chloride_to_carboxylic_acid_rxn:
            print("- No acyl chloride to carboxylic acid reaction")
        if not has_ester_throughout:
            print("- Ester group not preserved in key intermediates")
        if not correct_order:
            print("- Transformations not in correct order")

    return correct_sequence
