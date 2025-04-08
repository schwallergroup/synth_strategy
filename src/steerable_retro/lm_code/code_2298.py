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
    Detects if the route contains a late-stage esterification (carboxylic acid to ester)
    occurring in the final steps of the synthesis.
    """
    found_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal found_esterification

        if node["type"] == "reaction" and depth <= 2:  # Late stage = low depth
            try:
                if "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Checking reaction at depth {depth}: {rsmi}")

                    # Check if this is a hydrolysis reaction (esterification in retrosynthesis)
                    is_hydrolysis = checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )

                    if is_hydrolysis:
                        print(
                            f"Found hydrolysis reaction (esterification in retrosynthesis) at depth {depth}"
                        )
                        found_esterification = True
                        return

                    # Backup check: look for ester in reactants and carboxylic acid in product
                    # This corresponds to esterification in retrosynthesis
                    product_has_acid = checker.check_fg("Carboxylic acid", product)
                    reactant_has_ester = False
                    for r in reactants:
                        if checker.check_fg("Ester", r):
                            print(f"Found ester in reactant: {r}")
                            reactant_has_ester = True
                            break

                    if reactant_has_ester and product_has_acid:
                        print(
                            f"Found ester in reactant and carboxylic acid in product at depth {depth}"
                        )
                        found_esterification = True
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_esterification
