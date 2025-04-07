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
    This function detects if the synthetic route involves a Suzuki coupling reaction.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            # Extract reaction SMILES
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction: {rsmi}")

                # Check for Suzuki coupling using the checker function
                # Check all variants of Suzuki coupling
                suzuki_variants = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with boronic esters",
                ]

                for variant in suzuki_variants:
                    if checker.check_reaction(variant, rsmi):
                        print(f"Detected {variant}")
                        suzuki_detected = True
                        return

                # If the direct check didn't work, try a more detailed analysis
                if not suzuki_detected:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check for boronic acid or ester in reactants
                    has_boronic = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                            "Boronic ester", reactant
                        ):
                            has_boronic = True
                            print(f"Found boronic acid/ester in reactant: {reactant}")
                            break

                    # Check for aryl halide or triflate in reactants
                    has_aryl_halide_or_triflate = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Aromatic halide", reactant) or checker.check_fg(
                            "Triflate", reactant
                        ):
                            has_aryl_halide_or_triflate = True
                            print(f"Found aryl halide/triflate in reactant: {reactant}")
                            break

                    # Check for biaryl formation
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")
                    has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)

                    if has_boronic and has_aryl_halide_or_triflate and has_biaryl:
                        print("Detected Suzuki coupling based on reactant/product analysis")
                        suzuki_detected = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            if not suzuki_detected:  # Stop traversal if already detected
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return suzuki_detected
