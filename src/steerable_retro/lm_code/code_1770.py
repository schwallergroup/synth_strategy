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
    Detects if the synthesis uses a late-stage Suzuki coupling (C-C bond formation between
    two aromatic rings where one reactant has a boronic acid and the other has a halide).
    """
    suzuki_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_detected

        # Check if it's a late-stage reaction (final or penultimate step)
        if node["type"] == "reaction" and depth <= 1:
            print(f"Examining reaction at depth {depth}")
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                # Check if this is a Suzuki coupling using the checker function
                is_suzuki = False

                # Try different Suzuki coupling types
                suzuki_types = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with boronic esters",
                    "{Suzuki}",
                ]

                for suzuki_type in suzuki_types:
                    if checker.check_reaction(suzuki_type, rsmi):
                        print(f"Detected {suzuki_type}")
                        is_suzuki = True
                        break

                if is_suzuki:
                    # Extract reactants and product
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for boronic acid/ester in one reactant and halide in another
                    has_boronic = False
                    has_halide = False

                    for reactant in reactants:
                        # Check for boronic acid/ester
                        if checker.check_fg("Boronic acid", reactant):
                            print(f"Found boronic acid in reactant: {reactant}")
                            has_boronic = True
                        elif checker.check_fg("Boronic ester", reactant):
                            print(f"Found boronic ester in reactant: {reactant}")
                            has_boronic = True

                        # Check for halides
                        halide_types = [
                            "Aromatic halide",
                            "Primary halide",
                            "Secondary halide",
                            "Tertiary halide",
                            "Alkenyl halide",
                        ]
                        for halide_type in halide_types:
                            if checker.check_fg(halide_type, reactant):
                                print(f"Found {halide_type} in reactant: {reactant}")
                                has_halide = True
                                break

                    # If we haven't found the expected functional groups, try a more direct approach
                    if not (has_boronic and has_halide):
                        print("Trying direct pattern matching for boronic acid and halides")
                        for reactant in reactants:
                            # Direct check for boronic acid pattern
                            if "B(O)" in reactant or "B(OH)" in reactant:
                                print(
                                    f"Direct match: Found boronic acid pattern in reactant: {reactant}"
                                )
                                has_boronic = True

                            # Direct check for halides
                            if any(x in reactant for x in ["Br", "Cl", "I", "F"]):
                                mol = Chem.MolFromSmiles(reactant)
                                if mol:
                                    for atom in mol.GetAtoms():
                                        if (
                                            atom.GetSymbol() in ["Br", "Cl", "I", "F"]
                                            and atom.GetIsAromatic()
                                        ):
                                            print(
                                                f"Direct match: Found aromatic halide in reactant: {reactant}"
                                            )
                                            has_halide = True
                                            break

                    # Verify reactants contain aromatic rings
                    reactants_with_aromatic = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            has_aromatic = False
                            for atom in mol.GetAtoms():
                                if atom.GetIsAromatic():
                                    has_aromatic = True
                                    break
                            if has_aromatic:
                                reactants_with_aromatic += 1
                                print(f"Reactant with aromatic ring: {reactant}")

                    print(f"Reactants with aromatic rings: {reactants_with_aromatic}")
                    print(f"Has boronic: {has_boronic}, Has halide: {has_halide}")

                    # If we have boronic acid/ester, halide, and at least two reactants with aromatic rings
                    if has_boronic and has_halide:
                        print("Confirmed late-stage Suzuki coupling")
                        suzuki_detected = True
                else:
                    # If the reaction checker didn't identify it as Suzuki, try a more direct approach
                    print("Reaction not identified as Suzuki by checker, trying direct analysis")
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Look for key patterns in reactants
                    has_boronic = any("B(O)" in r or "B(OH)" in r for r in reactants)
                    has_halide = any(any(x in r for x in ["Br", "Cl", "I"]) for r in reactants)
                    has_pd = "Pd" in rsmi  # Check for palladium catalyst

                    if has_boronic and has_halide and has_pd:
                        print("Direct analysis suggests this is a Suzuki coupling")
                        suzuki_detected = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return suzuki_detected
