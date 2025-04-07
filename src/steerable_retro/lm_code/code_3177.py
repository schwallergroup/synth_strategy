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
    This function detects if the synthetic route preserves stereochemistry while
    modifying functional groups around a stereocenter.
    """
    # Track if we found preserved stereochemistry with functional group changes
    preserved_stereo = False

    # Track molecules at each depth to analyze changes
    molecules_by_depth = {}

    def dfs_traverse(node, depth=0):
        nonlocal preserved_stereo

        if node["type"] == "reaction":
            # Get reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                # Continue traversal even if no reaction SMILES
                for child in node.get("children", []):
                    dfs_traverse(child, depth + 1)
                return

            # Extract product
            product = rsmi.split(">")[-1]

            # Store molecule at this depth
            try:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    molecules_by_depth[depth] = {"mol": product_mol, "smiles": product}
            except Exception as e:
                print(f"Error creating RDKit molecule at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have molecules at multiple depths
    if len(molecules_by_depth) > 1:
        # Sort depths from early to late stage (high to low)
        depths = sorted(molecules_by_depth.keys(), reverse=True)

        # Check for stereochemistry preservation across transformations
        has_stereo = False

        # List of functional groups to check
        functional_groups = [
            "Nitro group",
            "Primary amine",
            "Secondary amine",
            "Tertiary amine",
            "Nitrile",
            "Carboxylic acid",
            "Ester",
            "Amide",
            "Alcohol",
            "Aldehyde",
            "Ketone",
            "Alkene",
            "Alkyne",
            "Halide",
        ]

        # Check if initial molecule has stereochemistry
        initial_depth = depths[0]  # Earliest stage
        initial_mol = molecules_by_depth[initial_depth]["mol"]
        initial_smiles = molecules_by_depth[initial_depth]["smiles"]

        stereo_centers = Chem.FindMolChiralCenters(initial_mol)
        if stereo_centers:
            has_stereo = True
            print(f"Found {len(stereo_centers)} stereocenters in initial molecule")

            # Track functional groups at each stereocenter
            for i in range(len(depths) - 1):
                early_depth = depths[i]
                late_depth = depths[i + 1]

                early_mol = molecules_by_depth[early_depth]["mol"]
                late_mol = molecules_by_depth[late_depth]["mol"]

                early_smiles = molecules_by_depth[early_depth]["smiles"]
                late_smiles = molecules_by_depth[late_depth]["smiles"]

                # Get stereocenters for both molecules
                early_stereo = dict(Chem.FindMolChiralCenters(early_mol))
                late_stereo = dict(Chem.FindMolChiralCenters(late_mol))

                # Check if stereocenters are preserved
                if len(early_stereo) > 0 and len(late_stereo) > 0:
                    # Check for functional group changes
                    fg_changes = False

                    for fg in functional_groups:
                        early_has_fg = checker.check_fg(fg, early_smiles)
                        late_has_fg = checker.check_fg(fg, late_smiles)

                        if early_has_fg != late_has_fg:
                            print(
                                f"Functional group change detected: {fg} - from {early_has_fg} to {late_has_fg}"
                            )
                            fg_changes = True

                    # If we have functional group changes and stereocenters in both molecules
                    if fg_changes and len(early_stereo) > 0 and len(late_stereo) > 0:
                        # Check if the reaction between these molecules preserves stereochemistry
                        # Get the reaction SMILES
                        reaction_found = False
                        for node in route.get("children", []):
                            if node["type"] == "reaction":
                                rsmi = node.get("metadata", {}).get("rsmi", "")
                                if rsmi:
                                    reactants = rsmi.split(">")[0].split(".")
                                    product = rsmi.split(">")[-1]

                                    # Check if this reaction connects our two molecules
                                    if (
                                        any(
                                            Chem.MolFromSmiles(r).GetNumAtoms()
                                            == early_mol.GetNumAtoms()
                                            for r in reactants
                                        )
                                        and Chem.MolFromSmiles(product).GetNumAtoms()
                                        == late_mol.GetNumAtoms()
                                    ):
                                        reaction_found = True

                                        # Check if stereocenters are preserved in this reaction
                                        if len(late_stereo) >= len(early_stereo):
                                            print(
                                                "Stereochemistry preserved with functional group changes"
                                            )
                                            preserved_stereo = True
                                            break

                        if not reaction_found:
                            # If we didn't find a direct reaction, just check if stereocenters are preserved
                            if len(late_stereo) >= len(early_stereo):
                                print("Stereochemistry preserved with functional group changes")
                                preserved_stereo = True

    return preserved_stereo
