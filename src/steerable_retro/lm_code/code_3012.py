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
    This function detects Suzuki coupling for biaryl formation.
    Looks for reactions where a boronic acid and aryl halide form a biaryl system.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # First, directly check if this is a Suzuki coupling reaction
            is_suzuki = False
            suzuki_reaction_types = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic esters",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with sulfonic esters",
                "{Suzuki}",  # General Suzuki reaction type
            ]

            for rxn_type in suzuki_reaction_types:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Detected {rxn_type}")
                    is_suzuki = True
                    break

            # If not detected by reaction type, check for characteristic patterns
            if not is_suzuki:
                # Check if we have a boronic compound and an aryl halide/triflate
                has_boronic_compound = False
                has_aryl_electrophile = False

                for reactant in reactants:
                    if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                        "Boronic ester", reactant
                    ):
                        print(f"Found boronic compound in reactant: {reactant}")
                        has_boronic_compound = True

                    if checker.check_fg("Aromatic halide", reactant) or checker.check_fg(
                        "Triflate", reactant
                    ):
                        print(f"Found aryl electrophile in reactant: {reactant}")
                        has_aryl_electrophile = True

                # If we have both required components, this might be a Suzuki coupling
                if has_boronic_compound and has_aryl_electrophile:
                    print(
                        "Found both boronic compound and aryl electrophile, checking product for biaryl system"
                    )

                    # Check if the product has a biaryl system
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Better pattern for biaryl detection
                        biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")
                        if product_mol.HasSubstructMatch(biaryl_pattern):
                            print(
                                f"Detected potential Suzuki coupling with biaryl formation: {rsmi}"
                            )
                            is_suzuki = True

            if is_suzuki:
                # Verify we have the right reactants
                has_boronic_compound = False
                has_aryl_electrophile = False

                for reactant in reactants:
                    print(f"Checking reactant: {reactant}")
                    # Check for boronic acids or esters
                    if checker.check_fg("Boronic acid", reactant):
                        print("Found boronic acid")
                        has_boronic_compound = True
                    elif checker.check_fg("Boronic ester", reactant):
                        print("Found boronic ester")
                        has_boronic_compound = True

                    # Check for aryl halides or triflates
                    if checker.check_fg("Aromatic halide", reactant):
                        print("Found aromatic halide")
                        has_aryl_electrophile = True
                    elif checker.check_fg("Triflate", reactant):
                        print("Found triflate")
                        has_aryl_electrophile = True

                print(f"Has boronic compound: {has_boronic_compound}")
                print(f"Has aryl electrophile: {has_aryl_electrophile}")

                # Verify the product has a biaryl system
                if has_boronic_compound and has_aryl_electrophile:
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol:
                        # Better pattern for biaryl detection
                        biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")

                        if product_mol.HasSubstructMatch(biaryl_pattern):
                            print(f"Detected biaryl system in product: {product}")
                            suzuki_detected = True
                            print(f"Confirmed Suzuki coupling for biaryl formation: {rsmi}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected
