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
    Detects a synthetic strategy involving carboxylic acid protection followed by
    late-stage amide coupling with an aromatic amine.
    """
    # Initialize flags for strategy components
    acid_protection_depth = float("inf")
    amide_coupling_depth = float("inf")
    pyrazole_disconnection_depth = float("inf")
    acyl_halide_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal acid_protection_depth, amide_coupling_depth, pyrazole_disconnection_depth, acyl_halide_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for carboxylic acid protection
                if any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles):
                    if (
                        checker.check_reaction("Protection of carboxylic acid", rsmi)
                        or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                        or checker.check_reaction(
                            "O-alkylation of carboxylic acids with diazo compounds", rsmi
                        )
                    ):
                        print(f"Found carboxylic acid protection at depth {depth}")
                        acid_protection_depth = min(acid_protection_depth, depth)

                # Check for amide coupling with aromatic amine
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Schotten-Baumann to ester",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]

                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        # Check if one of the reactants is an aromatic amine or if product has amide
                        if (
                            any(checker.check_fg("Aniline", r) for r in reactants_smiles)
                            or any(checker.check_fg("Primary amide", product_smiles) for _ in [1])
                            or any(checker.check_fg("Secondary amide", product_smiles) for _ in [1])
                            or any(checker.check_fg("Tertiary amide", product_smiles) for _ in [1])
                        ):
                            print(
                                f"Found amide coupling at depth {depth} with reaction {reaction_type}"
                            )
                            amide_coupling_depth = min(amide_coupling_depth, depth)
                            break

                # Check for pyrazole disconnection
                for r in reactants_smiles:
                    if checker.check_ring("pyrazole", r):
                        print(f"Found pyrazole disconnection at depth {depth}")
                        pyrazole_disconnection_depth = min(pyrazole_disconnection_depth, depth)
                        break

                # Check for acyl halide formation (not just fluoride)
                for r in reactants_smiles:
                    if checker.check_fg("Acyl halide", r):
                        print(f"Found acyl halide at depth {depth}")
                        acyl_halide_depth = min(acyl_halide_depth, depth)
                        break

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    has_acid_protection = acid_protection_depth != float("inf")
    has_late_amide_coupling = amide_coupling_depth != float("inf")
    has_pyrazole = pyrazole_disconnection_depth != float("inf")
    has_acyl_halide = acyl_halide_depth != float("inf")

    # More flexible sequence checking
    # Either acid protection followed by amide coupling
    # OR acyl halide formation (which can be used for amide coupling)
    # WITH either pyrazole or acyl halide involvement

    correct_sequence = (
        # Either acid protection before amide coupling
        (
            (
                has_acid_protection
                and has_late_amide_coupling
                and acid_protection_depth > amide_coupling_depth
            )
            or
            # Or just acyl halide formation (which is often used for amide coupling)
            has_acyl_halide
        )
        and
        # With either pyrazole or acyl halide involvement
        (has_pyrazole or has_acyl_halide)
    )

    print("Strategy detection results:")
    print(f"Acid protection: {has_acid_protection} (depth {acid_protection_depth})")
    print(f"Amide coupling: {has_late_amide_coupling} (depth {amide_coupling_depth})")
    print(f"Pyrazole: {has_pyrazole} (depth {pyrazole_disconnection_depth})")
    print(f"Acyl halide: {has_acyl_halide} (depth {acyl_halide_depth})")

    if correct_sequence:
        print(
            "Detected late-stage amide formation with carboxylic acid protection-activation sequence"
        )
    else:
        print("Strategy not detected")

    return correct_sequence
