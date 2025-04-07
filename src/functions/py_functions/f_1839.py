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
    Detects if the synthesis route includes Friedel-Crafts acylation.
    Looks for addition of an acyl group to an aromatic ring.
    """
    acylation_found = False

    def dfs_traverse(node):
        nonlocal acylation_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Use the checker function to directly identify Friedel-Crafts acylation
            if checker.check_reaction("Friedel-Crafts acylation", rsmi):
                print(f"Friedel-Crafts acylation detected: {rsmi}")
                acylation_found = True

            # As a fallback, check for the reaction components manually
            if not acylation_found:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has a ketone attached to an aromatic ring
                if checker.check_fg("Ketone", product) and any(
                    checker.check_ring(ring, product)
                    for ring in ["benzene", "naphthalene", "anthracene"]
                ):

                    # Check if one reactant is an acyl halide and another is an aromatic compound
                    has_acyl_halide = any(
                        checker.check_fg("Acyl halide", reactant)
                        for reactant in reactants
                    )
                    has_aromatic = any(
                        any(
                            checker.check_ring(ring, reactant)
                            for ring in ["benzene", "naphthalene", "anthracene"]
                        )
                        for reactant in reactants
                    )

                    if has_acyl_halide and has_aromatic:
                        print(f"Friedel-Crafts acylation components detected: {rsmi}")
                        acylation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return acylation_found
