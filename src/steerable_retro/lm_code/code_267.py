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
    This function detects a strategy involving sequential nitrogen functional group
    transformations (nitro → amine → hydrazine → heterocycle).
    """
    # Track transformation sequence with molecule SMILES for reference
    transformation_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for functional groups in molecules
            if checker.check_fg("Nitro group", mol_smiles):
                transformation_sequence.append(("nitro", depth, mol_smiles))

            if checker.check_fg("Primary amine", mol_smiles):
                transformation_sequence.append(("amine", depth, mol_smiles))

            if checker.check_fg("Hydrazine", mol_smiles):
                transformation_sequence.append(("hydrazine", depth, mol_smiles))

            if checker.check_ring("pyrazole", mol_smiles):
                transformation_sequence.append(("pyrazole", depth, mol_smiles))

        elif node["type"] == "reaction":
            # Check for specific reaction types
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction to amine
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction reaction at depth {depth}")
                    transformation_sequence.append(("nitro_to_amine", depth, rsmi))
                # Alternative check for nitro reduction
                elif (
                    any(checker.check_fg("Nitro group", r) for r in reactants)
                    and checker.check_fg("Primary amine", product)
                    and not any(checker.check_fg("Nitro group", product))
                ):
                    print(f"Found nitro reduction (by FG check) at depth {depth}")
                    transformation_sequence.append(("nitro_to_amine", depth, rsmi))

                # Check for amine to hydrazine conversion
                if any(
                    checker.check_fg("Primary amine", r) for r in reactants
                ) and checker.check_fg("Hydrazine", product):
                    print(f"Found amine to hydrazine conversion at depth {depth}")
                    transformation_sequence.append(("amine_to_hydrazine", depth, rsmi))

                # Check for hydrazine to pyrazole formation
                if any(checker.check_fg("Hydrazine", r) for r in reactants) and checker.check_ring(
                    "pyrazole", product
                ):
                    print(f"Found hydrazine to pyrazole formation at depth {depth}")
                    transformation_sequence.append(("hydrazine_to_pyrazole", depth, rsmi))

                # Alternative check for pyrazole formation
                elif checker.check_reaction("Pyrazole formation", rsmi):
                    print(f"Found pyrazole formation reaction at depth {depth}")
                    transformation_sequence.append(("hydrazine_to_pyrazole", depth, rsmi))

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth (descending) to get chronological order
    transformation_sequence.sort(key=lambda x: x[1], reverse=True)

    # Extract just the transformation types
    sequence = [item[0] for item in transformation_sequence]
    print(f"Transformation sequence: {sequence}")

    # Check for the presence of key transformations
    has_nitro = "nitro" in sequence
    has_amine = "amine" in sequence
    has_hydrazine = "hydrazine" in sequence
    has_pyrazole = "pyrazole" in sequence

    # Check for specific transformation reactions
    has_nitro_to_amine = "nitro_to_amine" in sequence
    has_amine_to_hydrazine = "amine_to_hydrazine" in sequence
    has_hydrazine_to_pyrazole = "hydrazine_to_pyrazole" in sequence

    # Check for the correct sequence order
    correct_order = False

    # Get indices of each functional group
    nitro_indices = [i for i, x in enumerate(sequence) if x == "nitro"]
    amine_indices = [i for i, x in enumerate(sequence) if x == "amine"]
    hydrazine_indices = [i for i, x in enumerate(sequence) if x == "hydrazine"]
    pyrazole_indices = [i for i, x in enumerate(sequence) if x == "pyrazole"]

    # Check if we have all functional groups and they appear in the correct order
    if nitro_indices and amine_indices and hydrazine_indices and pyrazole_indices:
        # Check if there exists any valid sequence where nitro -> amine -> hydrazine -> pyrazole
        for ni in nitro_indices:
            for ai in amine_indices:
                if ai > ni:  # amine comes after nitro
                    for hi in hydrazine_indices:
                        if hi > ai:  # hydrazine comes after amine
                            for pi in pyrazole_indices:
                                if pi > hi:  # pyrazole comes after hydrazine
                                    correct_order = True
                                    print(
                                        f"Correct sequence order detected: nitro({ni}) -> amine({ai}) -> hydrazine({hi}) -> pyrazole({pi})"
                                    )
                                    break
                            if correct_order:
                                break
                    if correct_order:
                        break
            if correct_order:
                break

    # Check for reaction transformations
    reaction_sequence = has_nitro_to_amine and has_amine_to_hydrazine and has_hydrazine_to_pyrazole

    # Check for functional group sequence
    fg_sequence = has_nitro and has_amine and has_hydrazine and has_pyrazole and correct_order

    # Determine if the strategy is present
    # We need either specific transformation reactions OR the correct sequence of functional groups
    strategy_present = reaction_sequence or fg_sequence

    print(f"Sequential nitrogen transformations strategy detected: {strategy_present}")
    print(f"Reaction sequence: {reaction_sequence}, FG sequence: {fg_sequence}")

    return strategy_present
