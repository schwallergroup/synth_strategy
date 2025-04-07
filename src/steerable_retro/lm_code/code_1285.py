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
    This function detects phosphonate formation via Arbuzov reaction.
    It looks for a reaction where a benzyl halide reacts with a phosphite to form a phosphonate.
    """
    arbuzov_reaction_found = False

    def dfs_traverse(node):
        nonlocal arbuzov_reaction_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains phosphonate pattern
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Phosphonate pattern: P(=O)(OR)₂R
                    phosphonate_pattern = Chem.MolFromSmarts("[P](=[O])([O,C])([O,C])[C,c]")
                    if product_mol.HasSubstructMatch(phosphonate_pattern):
                        print(f"Product contains phosphonate group: {product}")

                        # Check for benzyl/alkyl halide in reactants
                        has_halide = False
                        has_phosphite = False

                        for reactant in reactants:
                            # Check for primary, secondary, or tertiary halides
                            if (
                                checker.check_fg("Primary halide", reactant)
                                or checker.check_fg("Secondary halide", reactant)
                                or checker.check_fg("Tertiary halide", reactant)
                            ):
                                has_halide = True
                                print(f"Found halide in reactants: {reactant}")

                            # Check for phosphite structure
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # Phosphite pattern: P(OR)₃
                                phosphite_pattern = Chem.MolFromSmarts(
                                    "[P]([O][C,c])([O][C,c])([O][C,c])"
                                )
                                if reactant_mol.HasSubstructMatch(phosphite_pattern):
                                    has_phosphite = True
                                    print(f"Found phosphite in reactants: {reactant}")

                        # Verify Arbuzov reaction: halide + phosphite → phosphonate
                        if has_halide and has_phosphite:
                            print(
                                "Arbuzov reaction detected: halide and phosphite forming phosphonate"
                            )
                            arbuzov_reaction_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return arbuzov_reaction_found
