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
    This function detects if the synthesis route involves incorporation of a terminal alkene-containing fragment.
    Terminal alkenes include vinyl (C=CH2) and allyl (CH2=CH-CH2-) groups.
    """
    terminal_alkene_incorporated = False

    def dfs_traverse(node, depth=0):
        nonlocal terminal_alkene_incorporated

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            # Check for terminal alkene in reactants
            terminal_alkene_reactant = None
            for r in reactants:
                if checker.check_fg("Vinyl", r) or checker.check_fg("Allyl", r):
                    print(f"Terminal alkene found in reactant at depth {depth}: {r}")
                    terminal_alkene_reactant = r
                    break

            # If terminal alkene found in reactants, check if it's incorporated
            if terminal_alkene_reactant:
                # Check if the reaction is a known alkene incorporation reaction
                if (
                    checker.check_reaction("Diels-Alder", rsmi)
                    or checker.check_reaction("Heck terminal vinyl", rsmi)
                    or checker.check_reaction("Heck_terminal_vinyl", rsmi)
                    or checker.check_reaction("thiol-ene reaction", rsmi)
                    or checker.check_reaction("Michael addition", rsmi)
                    or checker.check_reaction("aza-Michael addition primary", rsmi)
                    or checker.check_reaction("aza-Michael addition secondary", rsmi)
                    or checker.check_reaction("aza-Michael addition aromatic", rsmi)
                    or checker.check_reaction("thia-Michael addition", rsmi)
                    or checker.check_reaction("oxa-Michael addition", rsmi)
                ):
                    print(f"Terminal alkene incorporation reaction detected at depth {depth}")
                    terminal_alkene_incorporated = True
                else:
                    # Check if product no longer has the terminal alkene
                    # This indicates the alkene was incorporated/consumed
                    product_has_terminal_alkene = checker.check_fg(
                        "Vinyl", product_part
                    ) or checker.check_fg("Allyl", product_part)

                    # If there are multiple reactants and the product doesn't have a terminal alkene,
                    # it's likely the terminal alkene was incorporated
                    if len(reactants) > 1 and not product_has_terminal_alkene:
                        print(
                            f"Terminal alkene incorporation detected at depth {depth} - alkene consumed"
                        )
                        terminal_alkene_incorporated = True
                    # If there's only one reactant with terminal alkene and the product structure changed
                    elif len(reactants) == 1 and not product_has_terminal_alkene:
                        print(f"Terminal alkene transformation detected at depth {depth}")
                        terminal_alkene_incorporated = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return terminal_alkene_incorporated
