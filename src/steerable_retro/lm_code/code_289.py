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
    This function detects if the synthesis route includes C-Si bond formation,
    specifically aryl-Si bond formation from aryl halide.
    """
    has_c_si_formation = False

    def dfs_traverse(node):
        nonlocal has_c_si_formation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a known C-Si bond forming reaction
                if checker.check_reaction("Hiyama-Denmark Coupling", rsmi):
                    print(f"Detected C-Si bond formation via Hiyama-Denmark Coupling: {rsmi}")
                    has_c_si_formation = True
                    return

                # Check for alcohol protection with silyl ethers
                if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                    print(f"Detected C-Si bond formation via silyl protection: {rsmi}")
                    has_c_si_formation = True
                    return

                # Check for other C-Si bond forming reactions
                # Check if reactants contain aryl halide
                has_aryl_halide = any(
                    checker.check_fg("Aromatic halide", reactant) for reactant in reactants
                )

                # Check for silicon reagents in reactants (more comprehensive check)
                has_si_reagent = any("[Si" in reactant for reactant in reactants)

                # Check if product contains silyl groups
                has_silyl_group = checker.check_fg(
                    "Silyl protective group", product
                ) or checker.check_fg("Silane", product)

                # Check if product contains C-Si bond
                product_mol = Chem.MolFromSmiles(product)
                has_c_si_bond = False
                if product_mol is not None:
                    # Try to find C-Si bond in product
                    if product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("c-[Si]")
                    ) or product_mol.HasSubstructMatch(Chem.MolFromSmarts("C-[Si]")):
                        has_c_si_bond = True

                # Check for aryl-Si formation from aryl halide
                if has_aryl_halide and has_si_reagent and (has_silyl_group or has_c_si_bond):
                    print(f"Detected C-Si bond formation: {rsmi}")
                    has_c_si_formation = True

                # Check for TMS ether formation
                if (
                    any(
                        checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                        for reactant in reactants
                    )
                    and has_si_reagent
                    and checker.check_fg("TMS ether protective group", product)
                ):
                    print(f"Detected C-Si bond formation via TMS protection: {rsmi}")
                    has_c_si_formation = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_c_si_formation
