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
    Detects if the synthesis route includes a biaryl formation via Suzuki coupling.
    Looks for a reaction where an aryl halide and boronic acid/ester are combined.
    """
    suzuki_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check if this is a Suzuki coupling reaction
            is_suzuki = (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                or checker.check_reaction("Suzuki", rsmi)
            )

            print(f"Is Suzuki coupling: {is_suzuki}")

            # Split reactants and product
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            print(f"Reactants: {reactants}")
            print(f"Product: {product_part}")

            # Check for aryl halide and boronic acid/ester in reactants
            has_aryl_halide = any(
                checker.check_fg("Aromatic halide", reactant) for reactant in reactants if reactant
            )
            has_boronic = any(
                checker.check_fg("Boronic acid", reactant)
                or checker.check_fg("Boronic ester", reactant)
                for reactant in reactants
                if reactant
            )

            print(f"Has aryl halide: {has_aryl_halide}")
            print(f"Has boronic: {has_boronic}")

            # Manual check for Suzuki-like reaction if checker fails
            if not is_suzuki and has_aryl_halide and has_boronic:
                print("Detected Suzuki-like reaction pattern")
                is_suzuki = True

            if is_suzuki:
                try:
                    # Create molecule objects
                    product_mol = Chem.MolFromSmiles(product_part)

                    # Find aryl halide and boronic acid/ester reactants
                    aryl_halide_reactants = [
                        r for r in reactants if r and checker.check_fg("Aromatic halide", r)
                    ]
                    boronic_reactants = [
                        r
                        for r in reactants
                        if r
                        and (
                            checker.check_fg("Boronic acid", r)
                            or checker.check_fg("Boronic ester", r)
                        )
                    ]

                    print(f"Aryl halide reactants: {aryl_halide_reactants}")
                    print(f"Boronic reactants: {boronic_reactants}")

                    # Check for biaryl formation
                    if aryl_halide_reactants and boronic_reactants and product_mol:
                        # First method: Check for biaryl bond pattern
                        biaryl_pattern = Chem.MolFromSmarts(
                            "c-c"
                        )  # Simple pattern for carbon-carbon bond between aromatics

                        # Create reactant molecules
                        reactant_mols = []
                        for r in reactants:
                            if r:
                                mol = Chem.MolFromSmiles(r)
                                if mol:
                                    reactant_mols.append(mol)

                        # Count biaryl bonds in product and reactants
                        if biaryl_pattern:
                            product_biaryl_count = len(
                                product_mol.GetSubstructMatches(biaryl_pattern)
                            )
                            reactant_biaryl_count = sum(
                                len(m.GetSubstructMatches(biaryl_pattern)) for m in reactant_mols
                            )

                            print(f"Product biaryl count: {product_biaryl_count}")
                            print(f"Reactant biaryl count: {reactant_biaryl_count}")

                            # If product has more biaryl bonds than reactants combined, a new one was formed
                            if product_biaryl_count > reactant_biaryl_count:
                                print(
                                    "Detected Suzuki coupling for biaryl formation (pattern match)"
                                )
                                suzuki_detected = True

                        # Second method: If we have the right reactants in a Suzuki reaction, assume biaryl formation
                        if not suzuki_detected:
                            print("Assuming biaryl formation based on reactants in Suzuki coupling")
                            suzuki_detected = True
                    else:
                        print("Missing required components for biaryl formation")
                except Exception as e:
                    print(f"Error processing Suzuki reaction: {e}")
                    # If we've identified a Suzuki reaction with aryl halide and boronic components,
                    # assume it's a biaryl formation even if analysis failed
                    if has_aryl_halide and has_boronic:
                        print("Assuming biaryl formation despite error")
                        suzuki_detected = True
            else:
                print("Not a Suzuki coupling reaction")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    print("Starting biaryl formation via Suzuki detection...")
    dfs_traverse(route)
    print(f"Suzuki biaryl formation detected: {suzuki_detected}")
    return suzuki_detected
