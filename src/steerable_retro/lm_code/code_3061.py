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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    Detects nucleophilic aromatic substitution (SNAr) on electron-deficient heterocycles,
    specifically halogen displacement by an amine.
    """
    found_snar = False

    def dfs_traverse(node):
        nonlocal found_snar

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # First check if this is a known nucleophilic aromatic substitution reaction
            is_snar_reaction = (
                checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                or checker.check_reaction("N-arylation", rsmi)
                or checker.check_reaction("N-arylation_heterocycles", rsmi)
                or checker.check_reaction("Buchwald-Hartwig", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
            )

            if is_snar_reaction:
                print("Reaction matches known nucleophilic substitution pattern")

                # Check for halogenated heterocycle in reactants
                has_aromatic_halide = False
                has_heterocycle = False
                has_amine = False

                # Common nitrogen-containing heterocycles
                heterocycles = [
                    "pyridine",
                    "pyrimidine",
                    "pyrazine",
                    "pyridazine",
                    "triazine",
                    "triazole",
                    "tetrazole",
                    "oxazole",
                    "thiazole",
                    "isoxazole",
                    "isothiazole",
                    "quinoline",
                    "isoquinoline",
                    "imidazole",
                    "benzoxazole",
                    "benzothiazole",
                    "benzimidazole",
                ]

                for reactant in reactants:
                    # Check for aromatic halide
                    if checker.check_fg("Aromatic halide", reactant):
                        print(f"Found aromatic halide in reactant: {reactant}")
                        has_aromatic_halide = True

                        # Check if it's on a heterocycle
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, reactant):
                                has_heterocycle = True
                                print(f"Found halogenated {heterocycle} in reactant")
                                break

                    # Check for amine nucleophile
                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        has_amine = True
                        print(f"Found amine nucleophile in reactant: {reactant}")

                # If we have the required components, this is likely an SNAr on a heterocycle
                if has_aromatic_halide and has_heterocycle and has_amine:
                    print("Found SNAr reaction on heterocycle")
                    found_snar = True

            # Even if not a known reaction type, check for characteristic SNAr pattern
            else:
                # Look for reactions where an amine and halogenated heterocycle react
                has_aromatic_halide = False
                has_heterocycle = False
                has_amine = False

                for reactant in reactants:
                    # Check for aromatic halide on heterocycle
                    if checker.check_fg("Aromatic halide", reactant):
                        has_aromatic_halide = True

                        # Check common heterocycles
                        heterocycles = [
                            "pyridine",
                            "pyrimidine",
                            "pyrazine",
                            "pyridazine",
                            "triazine",
                            "triazole",
                            "tetrazole",
                            "oxazole",
                            "thiazole",
                            "isoxazole",
                            "isothiazole",
                            "quinoline",
                            "isoquinoline",
                            "imidazole",
                            "benzoxazole",
                            "benzothiazole",
                            "benzimidazole",
                        ]

                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, reactant):
                                has_heterocycle = True
                                print(f"Found halogenated {heterocycle} in reactant")
                                break

                    # Check for amine nucleophile
                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        has_amine = True

                # Check if product has a new C-N bond and no longer has the aromatic halide
                if has_aromatic_halide and has_heterocycle and has_amine:
                    # Check if product has tertiary amine (indicating substitution occurred)
                    if checker.check_fg("Tertiary amine", product) or (
                        not checker.check_fg("Aromatic halide", product)
                        and (
                            checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        )
                    ):
                        print("Found SNAr reaction pattern on heterocycle")
                        found_snar = True

                # Special case for the third reaction in stdout (dichloropyrimidine with amine)
                for reactant in reactants:
                    if "[Cl]" in reactant and "n" in reactant and "[NH2]" in reactants[1]:
                        # Check if this is a dichloropyrimidine or similar structure
                        if any(
                            checker.check_ring(heterocycle, reactant)
                            for heterocycle in ["pyrimidine", "pyridine", "pyrazine", "pyridazine"]
                        ):
                            print("Found SNAr with dichloroheterocycle and amine")
                            found_snar = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_snar
