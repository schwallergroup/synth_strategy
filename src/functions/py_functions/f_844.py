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
    Detects if the final product contains a trifluoromethyl group,
    which is introduced via a coupling reaction.
    """
    has_cf3_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_cf3_coupling

        if node["type"] == "mol" and depth == 0:  # Final product
            final_product_smiles = node["smiles"]
            print(f"Analyzing final product: {final_product_smiles}")

            # Check if final product contains CF3 group
            if not checker.check_fg("Trifluoro group", final_product_smiles):
                print("Final product does not contain CF3 group")
                return

            print("Final product contains CF3 group")

            # Check if any of the reactions leading to this product is a coupling reaction
            # that introduced the CF3 group
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    try:
                        rsmi = child["metadata"]["rsmi"]
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        print(f"Checking reaction: {rsmi}")

                        # Check if this is a coupling reaction
                        coupling_reactions = [
                            "Suzuki coupling with boronic acids",
                            "Suzuki coupling with boronic acids OTf",
                            "Suzuki coupling with boronic esters",
                            "Suzuki coupling with boronic esters OTf",
                            "Suzuki coupling with sulfonic esters",
                            "Negishi coupling",
                            "Stille reaction_aryl",
                            "Stille reaction_vinyl",
                            "Stille reaction_benzyl",
                            "Stille reaction_allyl",
                            "Stille reaction_aryl OTf",
                            "Stille reaction_vinyl OTf",
                            "Sonogashira alkyne_aryl halide",
                            "Sonogashira acetylene_aryl halide",
                            "Sonogashira alkyne_aryl OTf",
                            "Sonogashira acetylene_aryl OTf",
                            "Hiyama-Denmark Coupling",
                            "Kumada cross-coupling",
                            "Buchwald-Hartwig",
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                            "{Suzuki}",
                        ]

                        is_coupling = False
                        for rxn_type in coupling_reactions:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(f"Detected coupling reaction: {rxn_type}")
                                is_coupling = True
                                break

                        if not is_coupling:
                            # Check if it's a general coupling by looking for key patterns
                            if "Pd" in rsmi or "Ni" in rsmi:
                                print(
                                    "Detected potential coupling reaction with Pd or Ni catalyst"
                                )
                                is_coupling = True
                            else:
                                print("Not a coupling reaction")
                                continue

                        print("Found a coupling reaction")

                        # Check if any reactant contains CF3 group
                        for i, reactant in enumerate(reactants):
                            if checker.check_fg("Trifluoro group", reactant):
                                print(f"Found CF3 group in reactant {i}: {reactant}")
                                has_cf3_coupling = True
                                return

                        print("No reactant contains CF3 group")
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_cf3_coupling
