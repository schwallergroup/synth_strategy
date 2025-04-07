#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects pyrazole synthesis via enamine-hydrazine cyclization
    followed by sequential functionalization (halogenation and N-acylation).
    """
    # Initialize flags for key transformations
    has_enamine_formation = False
    has_pyrazole_formation = False
    has_halogenation = False
    has_n_acylation = False

    # Track the sequence of transformations (in retrosynthetic order)
    transformation_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal has_enamine_formation, has_pyrazole_formation, has_halogenation, has_n_acylation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Depth {depth}, Reaction: {rsmi}")

                # Check for pyrazole formation
                if (
                    checker.check_ring("pyrazole", product_smiles)
                    and any(checker.check_fg("Hydrazine", r) for r in reactants_smiles)
                    and not any(checker.check_ring("pyrazole", r) for r in reactants_smiles)
                ):
                    has_pyrazole_formation = True
                    transformation_sequence.append(("pyrazole_formation", depth))
                    print(f"Detected pyrazole formation at depth {depth}")

                # Check for enamine formation
                if (
                    checker.check_fg("Enamine", product_smiles)
                    and any(checker.check_fg("Ketone", r) for r in reactants_smiles)
                    and not any(checker.check_fg("Enamine", r) for r in reactants_smiles)
                ):
                    has_enamine_formation = True
                    transformation_sequence.append(("enamine_formation", depth))
                    print(f"Detected enamine formation at depth {depth}")

                # Check for halogenation of aromatic ring or pyrazole
                halogenation_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Fluorination",
                    "Chlorination",
                    "Bromination",
                    "Iodination",
                ]

                if any(checker.check_reaction(rxn, rsmi) for rxn in halogenation_reactions):
                    has_halogenation = True
                    transformation_sequence.append(("halogenation", depth))
                    print(f"Detected halogenation at depth {depth}")
                # Also check for halogenation by comparing reactants and products
                elif any(
                    checker.check_ring("pyrazole", r) for r in reactants_smiles
                ) and checker.check_ring("pyrazole", product_smiles):
                    # Check if a halogen is added to the pyrazole
                    reactant_pyrazole = next(
                        (r for r in reactants_smiles if checker.check_ring("pyrazole", r)), None
                    )
                    if reactant_pyrazole and product_smiles:
                        # Check if product has more halogens than reactant
                        if (
                            product_smiles.count("Cl") > reactant_pyrazole.count("Cl")
                            or product_smiles.count("Br") > reactant_pyrazole.count("Br")
                            or product_smiles.count("F") > reactant_pyrazole.count("F")
                            or product_smiles.count("I") > reactant_pyrazole.count("I")
                        ):
                            has_halogenation = True
                            transformation_sequence.append(("halogenation", depth))
                            print(f"Detected halogenation at depth {depth}")

                # Check for N-acylation
                acylation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                ]

                if any(checker.check_reaction(rxn, rsmi) for rxn in acylation_reactions):
                    has_n_acylation = True
                    transformation_sequence.append(("n_acylation", depth))
                    print(f"Detected N-acylation at depth {depth}")

        # Process children (retrosynthetic traversal)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth (retrosynthetic order)
    sorted_transformations = sorted(transformation_sequence, key=lambda x: x[1])
    transformation_types = [t[0] for t in sorted_transformations]

    print(f"Transformation sequence (retrosynthetic): {transformation_types}")

    # Check if the strategy is present
    # We need pyrazole formation, enamine formation, and either halogenation or N-acylation
    correct_sequence = False

    if has_pyrazole_formation and has_enamine_formation:
        if has_halogenation or has_n_acylation:
            correct_sequence = True

    print(f"Strategy detected: {correct_sequence}")
    print(f"Has pyrazole formation: {has_pyrazole_formation}")
    print(f"Has enamine formation: {has_enamine_formation}")
    print(f"Has halogenation: {has_halogenation}")
    print(f"Has N-acylation: {has_n_acylation}")

    return correct_sequence
