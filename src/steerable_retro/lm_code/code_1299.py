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
    This function detects if the synthesis involves an aromatic acylation step
    to form a diaryl ketone.
    """
    acylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal acylation_found

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a Friedel-Crafts acylation reaction
                if checker.check_reaction("Friedel-Crafts acylation", rsmi):
                    print(
                        f"Friedel-Crafts acylation reaction detected at depth {depth}, rsmi: {rsmi}"
                    )
                    acylation_found = True
                else:
                    # Alternative detection method
                    # Check if product contains a ketone and at least one aromatic ring
                    if checker.check_fg("Ketone", product):
                        # Check if product has aromatic rings
                        if checker.check_ring("benzene", product):
                            # Check if any reactant has acyl halide
                            acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                            # Check if any reactant has aromatic ring
                            aromatic = any(checker.check_ring("benzene", r) for r in reactants)

                            if acyl_halide and aromatic:
                                print(f"Aromatic acylation detected at depth {depth}, rsmi: {rsmi}")
                                print(f"Product: {product}, Reactants: {reactants}")
                                acylation_found = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return acylation_found
