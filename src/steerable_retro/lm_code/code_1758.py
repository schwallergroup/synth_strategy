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
    This function detects the transformation of a ketone to an oxime.
    """
    oxime_formation_found = False

    def dfs_traverse(node):
        nonlocal oxime_formation_found

        if node["type"] == "reaction":
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Analyzing reaction: {rsmi}")

                    # Check if this is a ketone to oxime transformation reaction
                    # First check if any reactant contains a ketone
                    ketone_in_reactants = False
                    for r in reactants:
                        try:
                            if checker.check_fg("Ketone", r):
                                ketone_in_reactants = True
                                print(f"Found ketone in reactant: {r}")
                                break
                        except Exception as e:
                            print(f"Error checking ketone in reactant {r}: {e}")

                    # Check if product contains an oxime
                    oxime_in_product = False
                    try:
                        if checker.check_fg("Oxime", product):
                            oxime_in_product = True
                            print(f"Found oxime in product: {product}")
                        else:
                            print(f"No oxime found in product: {product}")
                    except Exception as e:
                        print(f"Error checking oxime in product {product}: {e}")

                    # Check if this is a ketone to oxime transformation
                    if ketone_in_reactants and oxime_in_product:
                        # Check for hydroxylamine or its derivatives in reactants
                        hydroxylamine_present = False
                        for r in reactants:
                            try:
                                # More comprehensive check for hydroxylamine
                                if "[NH2][OH]" in r or "NH2OH" in r or r.upper() == "NH2OH":
                                    hydroxylamine_present = True
                                    print(f"Found hydroxylamine reagent: {r}")
                                    break
                                # Check for hydroxylamine pattern in atom-mapped SMILES
                                elif "[NH2:15][OH:16]" in r or "[NH2:15]O[H:16]" in r:
                                    hydroxylamine_present = True
                                    print(f"Found atom-mapped hydroxylamine reagent: {r}")
                                    break
                            except Exception as e:
                                print(f"Error checking hydroxylamine in reactant {r}: {e}")

                        # If we found hydroxylamine or the reaction involves ketone to oxime transformation
                        if hydroxylamine_present:
                            oxime_formation_found = True
                            print("Confirmed ketone to oxime transformation with hydroxylamine")
                        else:
                            # Check if any reactant contains NH2OH pattern
                            for r in reactants:
                                if "NH2" in r and "OH" in r:
                                    print(f"Found potential hydroxylamine pattern in: {r}")
                                    oxime_formation_found = True
                                    print(
                                        "Confirmed ketone to oxime transformation with potential hydroxylamine"
                                    )
                                    break
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: oxime_formation_found = {oxime_formation_found}")
    return oxime_formation_found
