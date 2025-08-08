#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects a synthetic strategy involving late-stage alkyne coupling
    between two complex fragments.
    """
    alkyne_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal alkyne_coupling_detected

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate reaction (late stage)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for specific alkyne coupling reactions
                sonogashira_coupling = (
                    checker.check_reaction("Sonogashira acetylene_aryl halide", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_aryl OTf", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl OTf", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_alkenyl halide", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_alkenyl halide", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_alkenyl OTf", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_alkenyl OTf", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_acyl halide", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_acyl halide", rsmi)
                )

                if sonogashira_coupling:
                    print(f"Sonogashira coupling detected at depth {depth}")
                    alkyne_coupling_detected = True
                    return

                # Check if product contains alkyne
                product_has_alkyne = checker.check_fg("Alkyne", product)

                if product_has_alkyne:
                    print(f"Product contains alkyne at depth {depth}")

                    # Check if we're joining two fragments with an alkyne bond
                    if len(reactants) >= 2:
                        # Look for alkyne in reactants
                        reactant_has_alkyne = [checker.check_fg("Alkyne", r) for r in reactants]

                        # Check for terminal alkyne specifically
                        terminal_alkyne_pattern = Chem.MolFromSmarts("C#C[H]")
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                        has_terminal_alkyne = [
                            mol and mol.HasSubstructMatch(terminal_alkyne_pattern)
                            for mol in reactant_mols
                            if mol
                        ]

                        # Check for aryl or vinyl halides
                        aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)
                        vinyl_halide = any(checker.check_fg("Alkenyl halide", r) for r in reactants)

                        print(f"Reactants with alkyne: {reactant_has_alkyne}")
                        print(f"Reactants with terminal alkyne: {has_terminal_alkyne}")
                        print(f"Aryl halide present: {aryl_halide}")
                        print(f"Vinyl halide present: {vinyl_halide}")

                        # If one reactant has terminal alkyne and another has aryl/vinyl halide,
                        # it's likely an alkyne coupling
                        if any(has_terminal_alkyne) and (aryl_halide or vinyl_halide):
                            print(f"Late-stage alkyne coupling detected at depth {depth}")
                            alkyne_coupling_detected = True
                            return

                        # If product has alkyne but no reactant does, it might be alkyne formation
                        if not any(reactant_has_alkyne) and product_has_alkyne:
                            print(f"Alkyne formation detected at depth {depth}")
                            alkyne_coupling_detected = True
                            return

                        # If one reactant has alkyne and another doesn't, check if it's a coupling
                        if any(reactant_has_alkyne) and not all(reactant_has_alkyne):
                            # Check if the reaction is likely a coupling
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol:
                                # Count carbon atoms in product and reactants
                                product_carbon_count = len(
                                    [
                                        atom
                                        for atom in product_mol.GetAtoms()
                                        if atom.GetAtomicNum() == 6
                                    ]
                                )
                                reactant_carbon_counts = [
                                    len(
                                        [
                                            atom
                                            for atom in mol.GetAtoms()
                                            if atom.GetAtomicNum() == 6
                                        ]
                                    )
                                    for mol in reactant_mols
                                    if mol
                                ]

                                # If product has approximately the sum of carbons from reactants,
                                # it's likely a coupling reaction
                                if product_carbon_count >= sum(reactant_carbon_counts) - 2:
                                    print(
                                        f"Late-stage alkyne coupling detected based on carbon count at depth {depth}"
                                    )
                                    alkyne_coupling_detected = True
                                    return

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return alkyne_coupling_detected
