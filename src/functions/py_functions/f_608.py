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
    This function detects direct C-H functionalization to form an aldehyde.
    """
    ch_to_aldehyde_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ch_to_aldehyde_detected

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains aldehyde
                if checker.check_fg("Aldehyde", product_smiles):
                    # Check if any reactant does NOT contain aldehyde
                    reactant_without_aldehyde = False
                    for reactant in reactants_smiles:
                        if not checker.check_fg("Aldehyde", reactant):
                            reactant_without_aldehyde = True
                            break

                    if reactant_without_aldehyde:
                        # Check for specific C-H functionalization reactions
                        if (
                            checker.check_reaction(
                                "Oxidation of alkene to aldehyde", rsmi
                            )
                            or checker.check_reaction(
                                "Oxidation of alcohol to aldehyde", rsmi
                            )
                            or checker.check_reaction("Aromatic hydroxylation", rsmi)
                            or checker.check_reaction(
                                "Directed ortho metalation of arenes", rsmi
                            )
                            or checker.check_reaction(
                                "Hydration of alkyne to aldehyde", rsmi
                            )
                            or checker.check_reaction(
                                "Bouveault aldehyde synthesis", rsmi
                            )
                            or checker.check_reaction(
                                "Homologation of aldehydes with formaldehyde", rsmi
                            )
                            or checker.check_reaction(
                                "Homologation of diazo compounds to aldehydes using formaldehyde",
                                rsmi,
                            )
                        ):

                            print(
                                f"Detected C-H functionalization to aldehyde at depth {depth}"
                            )
                            print(f"Reaction SMILES: {rsmi}")
                            ch_to_aldehyde_detected = True
                            return

                        # Check for primary alcohol oxidation (common C-H functionalization equivalent)
                        for reactant in reactants_smiles:
                            if checker.check_fg("Primary alcohol", reactant):
                                print(
                                    f"Detected primary alcohol oxidation to aldehyde at depth {depth}"
                                )
                                print(f"Reaction SMILES: {rsmi}")
                                ch_to_aldehyde_detected = True
                                return

                        # Additional check for general C-H functionalization
                        # by comparing atom counts and checking for aldehyde formation
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        for reactant in reactants_smiles:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and product_mol:
                                # Check if atom count difference is small (typical for C-H functionalization)
                                if (
                                    abs(
                                        reactant_mol.GetNumAtoms()
                                        - product_mol.GetNumAtoms()
                                    )
                                    <= 5
                                ):
                                    # Ensure we're not just detecting a simple functional group conversion
                                    if not (
                                        checker.check_fg("Ketone", reactant)
                                        or checker.check_fg("Carboxylic acid", reactant)
                                        or checker.check_fg("Ester", reactant)
                                        or checker.check_fg("Nitrile", reactant)
                                        or checker.check_fg("Amide", reactant)
                                        or checker.check_fg("Acyl halide", reactant)
                                    ):

                                        # Check for alkyne or alkene that could be oxidized to aldehyde
                                        if (
                                            checker.check_fg("Alkyne", reactant)
                                            or checker.check_fg("Alkene", reactant)
                                            or checker.check_fg("Allyl", reactant)
                                            or checker.check_fg("Vinyl", reactant)
                                        ):
                                            print(
                                                f"Detected potential C-H functionalization to aldehyde at depth {depth}"
                                            )
                                            print(f"Reaction SMILES: {rsmi}")
                                            ch_to_aldehyde_detected = True
                                            return
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return ch_to_aldehyde_detected
