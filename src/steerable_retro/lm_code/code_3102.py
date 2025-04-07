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
    Detects if the route includes an ester hydrolysis step (converting ester to carboxylic acid).
    """
    found_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal found_ester_hydrolysis

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an ester hydrolysis reaction
            if checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            ):
                # Verify that reactants contain ester and product contains carboxylic acid
                has_ester = any(checker.check_fg("Ester", reactant) for reactant in reactants)
                has_acid = checker.check_fg("Carboxylic acid", product)

                if has_ester and has_acid:
                    print(f"Found ester hydrolysis step: {rsmi}")
                    print(f"Reactants have ester: {has_ester}, Product has acid: {has_acid}")
                    found_ester_hydrolysis = True

            # Check for ester saponification reactions
            elif checker.check_reaction(
                "Ester saponification (methyl deprotection)", rsmi
            ) or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi):
                # Verify that reactants contain ester and product contains carboxylic acid
                has_ester = any(checker.check_fg("Ester", reactant) for reactant in reactants)
                has_acid = checker.check_fg("Carboxylic acid", product)

                if has_ester and has_acid:
                    print(f"Found ester saponification step: {rsmi}")
                    print(f"Reactants have ester: {has_ester}, Product has acid: {has_acid}")
                    found_ester_hydrolysis = True

            # Check for COOH ethyl deprotection (another form of ester hydrolysis)
            elif checker.check_reaction("COOH ethyl deprotection", rsmi):
                has_ester = any(checker.check_fg("Ester", reactant) for reactant in reactants)
                has_acid = checker.check_fg("Carboxylic acid", product)

                if has_ester and has_acid:
                    print(f"Found COOH ethyl deprotection step: {rsmi}")
                    found_ester_hydrolysis = True

            # Direct check for ester to carboxylic acid transformation
            else:
                has_ester = any(checker.check_fg("Ester", reactant) for reactant in reactants)
                has_acid = checker.check_fg("Carboxylic acid", product)
                has_acid_in_reactants = any(
                    checker.check_fg("Carboxylic acid", reactant) for reactant in reactants
                )

                # Make sure we're actually creating a new acid group, not just preserving an existing one
                if has_ester and has_acid and not has_acid_in_reactants:
                    print(f"Found ester to acid transformation: {rsmi}")
                    found_ester_hydrolysis = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_ester_hydrolysis
