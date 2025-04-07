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
    Detects a strategy involving sequential SNAr reactions on nitro-activated heterocycles.
    Looks for multiple reactions where an aromatic halide is replaced by a nitrogen nucleophile.
    """
    snar_reactions = 0
    snar_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal snar_reactions, snar_depths

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for SNAr reaction using reaction checkers
            is_snar = False

            # Check using known reaction patterns
            if (
                checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi)
                or checker.check_reaction("Goldberg coupling", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
            ):
                is_snar = True

            # If not identified by reaction checkers, check for characteristic patterns
            if not is_snar:
                # Check for aromatic halide in reactants
                has_aromatic_halide = any(
                    checker.check_fg("Aromatic halide", reactant) for reactant in reactants
                )

                # Check for amine nucleophile in reactants
                has_amine = any(
                    checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                    or checker.check_fg("Aniline", reactant)
                    for reactant in reactants
                )

                # Check for nitro activation in reactants
                has_nitro = any(checker.check_fg("Nitro group", reactant) for reactant in reactants)

                # Check for heterocycle in reactants
                has_heterocycle = any(
                    checker.check_ring(ring, reactant)
                    for ring in [
                        "pyridine",
                        "pyrimidine",
                        "pyrazine",
                        "pyridazine",
                        "triazine",
                        "furan",
                        "thiophene",
                        "pyrrole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                    ]
                    for reactant in reactants
                )

                # Check if product has new C-N bond
                prod_mol = Chem.MolFromSmiles(product)
                if (
                    has_aromatic_halide
                    and has_amine
                    and (has_nitro or has_heterocycle)
                    and prod_mol
                ):
                    # Check if the product has a new C-N bond where a halide was
                    if checker.check_fg("Aniline", product) or checker.check_fg(
                        "Secondary amine", product
                    ):
                        is_snar = True

            if is_snar:
                snar_reactions += 1
                snar_depths.append(depth)
                print(f"SNAr reaction detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 2 SNAr reactions
    if snar_reactions >= 2:
        # Check if they are sequential (within reasonable depth difference)
        snar_depths.sort()
        for i in range(len(snar_depths) - 1):
            if abs(snar_depths[i] - snar_depths[i + 1]) <= 2:
                print(
                    f"Sequential SNAr reactions detected at depths {snar_depths[i]} and {snar_depths[i+1]}"
                )
                return True

    return False
