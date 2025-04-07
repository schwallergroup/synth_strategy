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
    This function detects multiple amide couplings in the synthetic route.
    """
    amide_coupling_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check if this is an amide coupling reaction using the checker
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                ]

                is_amide_coupling = False
                for reaction_name in amide_coupling_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_amide_coupling = True
                        break

                # If not detected by reaction checker, try to detect by functional groups
                if not is_amide_coupling:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check for carboxylic acid and amine in reactants
                    has_carboxylic_acid = False
                    has_amine = False

                    for reactant_smiles in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", reactant_smiles):
                            has_carboxylic_acid = True

                        if checker.check_fg(
                            "Primary amine", reactant_smiles
                        ) or checker.check_fg("Secondary amine", reactant_smiles):
                            has_amine = True

                    # Check for amide in product
                    has_amide = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )

                    is_amide_coupling = has_carboxylic_acid and has_amine and has_amide

                if is_amide_coupling:
                    print(f"Amide coupling detected at depth {depth}")
                    amide_coupling_count += 1

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    print("Starting analysis for multiple amide coupling strategy...")
    dfs_traverse(route)
    print(f"Total amide couplings detected: {amide_coupling_count}")

    # Return True if at least 2 amide couplings are detected
    result = amide_coupling_count >= 2
    print(f"Multiple amide coupling strategy detected: {result}")
    return result
