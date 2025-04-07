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
    This function detects a synthetic strategy involving benzimidazole ring formation.
    """
    benzimidazole_formation_detected = False

    def has_o_phenylenediamine(smiles):
        """Check if molecule contains o-phenylenediamine structure"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Check for two adjacent aniline groups on benzene
        # First check if molecule has at least two aniline groups
        if smiles.count("N") < 2:
            return False

        # Check for ortho-diamine pattern on benzene
        pattern = Chem.MolFromSmarts("c1c(N)c(N)ccc1")
        if mol.HasSubstructMatch(pattern):
            return True

        # Also check for 1,2-diaminobenzene pattern
        pattern2 = Chem.MolFromSmarts("c1c(N)cc(N)cc1")
        return mol.HasSubstructMatch(pattern2)

    def dfs_traverse(node, depth=0):
        nonlocal benzimidazole_formation_detected

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains benzimidazole ring
                product_has_benzimidazole = checker.check_ring("benzimidazole", product_smiles)

                # Check if reactants contain benzimidazole ring
                reactants_have_benzimidazole = any(
                    checker.check_ring("benzimidazole", r) for r in reactants_smiles
                )

                # Only proceed if benzimidazole is in product but not in reactants (indicating formation)
                if product_has_benzimidazole and not reactants_have_benzimidazole:
                    print(f"Found potential benzimidazole formation at depth {depth}")

                    # Check for specific benzimidazole formation reactions from the provided list
                    if (
                        checker.check_reaction("benzimidazole formation from aldehyde", rsmi)
                        or checker.check_reaction("benzimidazole formation from acyl halide", rsmi)
                        or checker.check_reaction(
                            "benzimidazole formation from ester/carboxylic acid", rsmi
                        )
                        or checker.check_reaction(
                            "{benzimidazole_derivatives_carboxylic-acid/ester}", rsmi
                        )
                        or checker.check_reaction("{benzimidazole_derivatives_aldehyde}", rsmi)
                    ):
                        print(f"Detected benzimidazole formation reaction at depth {depth}")
                        benzimidazole_formation_detected = True
                        return

                    # Check for o-phenylenediamine in reactants
                    reactants_have_o_phenylenediamine = any(
                        has_o_phenylenediamine(r) for r in reactants_smiles
                    )

                    # Check for carbonyl compounds in reactants
                    reactants_have_aldehyde = any(
                        checker.check_fg("Aldehyde", r) for r in reactants_smiles
                    )
                    reactants_have_carboxylic = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )
                    reactants_have_ester = any(
                        checker.check_fg("Ester", r) for r in reactants_smiles
                    )
                    reactants_have_acyl_halide = any(
                        checker.check_fg("Acyl halide", r) for r in reactants_smiles
                    )
                    reactants_have_formaldehyde = any(
                        checker.check_fg("Formaldehyde", r) for r in reactants_smiles
                    )

                    # Check for combination of o-phenylenediamine and carbonyl compound
                    if reactants_have_o_phenylenediamine and (
                        reactants_have_aldehyde
                        or reactants_have_carboxylic
                        or reactants_have_ester
                        or reactants_have_acyl_halide
                        or reactants_have_formaldehyde
                    ):
                        print(
                            f"Detected benzimidazole formation from o-phenylenediamine and carbonyl at depth {depth}"
                        )
                        benzimidazole_formation_detected = True
                        return

                    # Additional check for two anilines in reactants
                    aniline_count = sum(
                        1 for r in reactants_smiles if checker.check_fg("Aniline", r)
                    )
                    if aniline_count >= 2 and (
                        reactants_have_aldehyde
                        or reactants_have_carboxylic
                        or reactants_have_ester
                        or reactants_have_acyl_halide
                        or reactants_have_formaldehyde
                    ):
                        print(
                            f"Detected potential benzimidazole formation from multiple anilines at depth {depth}"
                        )
                        benzimidazole_formation_detected = True
                        return
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Benzimidazole formation detected: {benzimidazole_formation_detected}")
    return benzimidazole_formation_detected
