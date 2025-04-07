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
    Detects if the synthesis route involves N-alkylation of a heterocycle.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is an N-arylation or N-alkylation reaction
                is_n_alkylation = (
                    checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction("N-arylation_heterocycles", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        rsmi,
                    )
                )

                # If not a known N-alkylation reaction, check for the pattern manually
                if not is_n_alkylation:
                    # Check for heterocycles in reactants
                    heterocycle_reactants = []
                    for r in reactants_smiles:
                        if (
                            checker.check_ring("pyrrole", r)
                            or checker.check_ring("imidazole", r)
                            or checker.check_ring("triazole", r)
                            or checker.check_ring("tetrazole", r)
                            or checker.check_ring("pyrazole", r)
                            or checker.check_ring("indole", r)
                            or checker.check_ring("benzimidazole", r)
                        ):
                            heterocycle_reactants.append(r)

                    # Check for alkyl halides in reactants
                    has_alkyl_halide = any(
                        checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        for r in reactants_smiles
                    )

                    # Check if product has N-alkylated heterocycle
                    if heterocycle_reactants and has_alkyl_halide:
                        # Check if the product contains an N-alkylated heterocycle
                        product_has_n_alkylated_heterocycle = (
                            checker.check_ring("pyrrole", product_smiles)
                            or checker.check_ring("imidazole", product_smiles)
                            or checker.check_ring("triazole", product_smiles)
                            or checker.check_ring("tetrazole", product_smiles)
                            or checker.check_ring("pyrazole", product_smiles)
                            or checker.check_ring("indole", product_smiles)
                            or checker.check_ring("benzimidazole", product_smiles)
                        )

                        if product_has_n_alkylated_heterocycle:
                            # Check for N-C bond in product that's not in reactants
                            product_mol = Chem.MolFromSmiles(product_smiles)
                            n_alkyl_pattern = Chem.MolFromSmarts("[n]-[C;!$(C=O)]")

                            if (
                                product_mol
                                and n_alkyl_pattern
                                and product_mol.HasSubstructMatch(n_alkyl_pattern)
                            ):
                                # Verify this pattern wasn't in all reactants
                                all_reactants_have_pattern = True
                                for r in heterocycle_reactants:
                                    r_mol = Chem.MolFromSmiles(r)
                                    if r_mol and not r_mol.HasSubstructMatch(
                                        n_alkyl_pattern
                                    ):
                                        all_reactants_have_pattern = False
                                        break

                                if not all_reactants_have_pattern:
                                    print(
                                        f"Detected N-alkylation of heterocycle in reaction: {rsmi}"
                                    )
                                    result = True
                else:
                    # If it's a known N-alkylation reaction, check if heterocycle is involved
                    has_heterocycle = any(
                        checker.check_ring("pyrrole", r)
                        or checker.check_ring("imidazole", r)
                        or checker.check_ring("triazole", r)
                        or checker.check_ring("tetrazole", r)
                        or checker.check_ring("pyrazole", r)
                        or checker.check_ring("indole", r)
                        or checker.check_ring("benzimidazole", r)
                        for r in reactants_smiles
                    )

                    if has_heterocycle:
                        print(
                            f"Detected N-alkylation of heterocycle via known reaction: {rsmi}"
                        )
                        result = True

            except Exception as e:
                print(f"Error in processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return result
