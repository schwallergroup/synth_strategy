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
    Detects if the synthesis uses a Grignard addition to a Weinreb amide.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is a Weinreb amide to ketone reaction
                if checker.check_reaction("Ketone from Weinreb amide", rsmi):
                    print(f"Found Ketone from Weinreb amide reaction at depth {depth}")

                    # Extract reactants to check for Grignard reagent
                    reactants = rsmi.split(">")[0].split(".")

                    # Check if any reactant is a Grignard reagent
                    for reactant in reactants:
                        if checker.check_fg("Magnesium halide", reactant):
                            print(f"Found Grignard reagent: {reactant}")
                            found_pattern = True
                            return

                # Alternative check: look for both components and verify reaction pattern
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                has_grignard = False
                has_weinreb = False

                for reactant in reactants:
                    # Check for Grignard reagent
                    if checker.check_fg("Magnesium halide", reactant):
                        print(f"Found Grignard reagent: {reactant}")
                        has_grignard = True

                    # Check for Weinreb amide structure
                    # Weinreb amides are tertiary amides with N-methoxy-N-methyl group
                    if checker.check_fg("Tertiary amide", reactant):
                        # More specific check for N-methoxy-N-methyl amide structure
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            for atom in mol.GetAtoms():
                                # Find nitrogen atoms
                                if atom.GetAtomicNum() == 7:
                                    # Check if connected to oxygen (methoxy group)
                                    for neighbor in atom.GetNeighbors():
                                        if neighbor.GetAtomicNum() == 8:
                                            print(f"Found potential Weinreb amide: {reactant}")
                                            has_weinreb = True

                # If both components are present and product is a ketone
                if has_grignard and has_weinreb:
                    if checker.check_fg("Ketone", product):
                        print(
                            f"Found Grignard-Weinreb pattern through component analysis at depth {depth}"
                        )
                        found_pattern = True
                    else:
                        print(f"Found Grignard and Weinreb components but product is not a ketone")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)
            if found_pattern:
                return

    # Start traversal from root
    dfs_traverse(route)
    return found_pattern
