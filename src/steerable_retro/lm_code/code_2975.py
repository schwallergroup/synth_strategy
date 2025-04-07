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
    Detects late-stage C-C bond formation via organometallic coupling (e.g., organozinc)
    in the final steps of the synthesis.
    """
    found_late_cc_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_cc_coupling

        # Print node type and depth for debugging
        print(f"Examining node type: {node['type']} at depth: {depth}")

        if node["type"] == "reaction" and depth <= 3:  # Late stage (expanded to depth 3)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction SMILES at depth {depth}: {rsmi}")

                # Check for specific C-C bond forming organometallic reactions
                organometallic_reactions = [
                    "Grignard_carbonyl",
                    "Grignard_alcohol",
                    "Negishi coupling",
                    "Stille reaction_aryl",
                    "Stille reaction_vinyl",
                    "Stille reaction_benzyl",
                    "Stille reaction_allyl",
                    "Stille reaction_vinyl OTf",
                    "Stille reaction_aryl OTf",
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with boronic esters OTf",
                    "Kumada cross-coupling",
                    "Hiyama-Denmark Coupling",
                    "Aryllithium cross-coupling",
                    "Grignard from aldehyde to alcohol",
                    "Grignard from ketone to alcohol",
                    "Grignard with CO2 to carboxylic acid",
                    "Grignard from nitrile to ketone",
                    "Olefination of ketones with Grignard reagents",
                    "Olefination of aldehydes with Grignard reagents",
                    "Preparation of trialkylsilanes with Grignard reagents",
                    "Formation of Grignard reagents",
                    "Reaction of alkyl halides with organometallic coumpounds",
                    "Preparation of organolithium compounds",
                    "Aerobic oxidation of Grignard reagents",
                    "Carboxylic acid from Li and CO2",
                    "Ketone from Li and CO2",
                    "Ketone from Li, Grignard and CO2",
                    "Ketone from Li, halide and CO2",
                    "Ketone from Grignard and CO2",
                    "Directed ortho metalation of arenes",
                ]

                for reaction_type in organometallic_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Found late-stage organometallic C-C coupling at depth {depth}: {reaction_type}"
                        )
                        found_late_cc_coupling = True
                        return  # Exit early once found

                # Check for organometallic functional groups in reactants
                if not found_late_cc_coupling:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check for organometallic functional groups in reactants
                        organometallic_fgs = [
                            "Magnesium halide",
                            "Zinc halide",
                            "Alkyl lithium",
                            "Aryl lithium",
                            "Tin",
                        ]
                        for reactant in reactants:
                            for fg in organometallic_fgs:
                                if checker.check_fg(fg, reactant):
                                    print(
                                        f"Found organometallic functional group {fg} in reactant at depth {depth}"
                                    )
                                    # Verify C-C bond formation by checking product
                                    found_late_cc_coupling = True
                                    return  # Exit early once found
                    except Exception as e:
                        print(f"Error processing reaction SMILES: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)
            if found_late_cc_coupling:
                return  # Exit early once found

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: {found_late_cc_coupling}")
    return found_late_cc_coupling
