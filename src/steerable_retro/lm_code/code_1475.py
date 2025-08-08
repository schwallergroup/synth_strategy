#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects a synthetic strategy where a nitrile group is preserved throughout the synthesis.
    """
    nitrile_present_in_final = False
    all_steps_preserve_nitrile = True

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_present_in_final, all_steps_preserve_nitrile

        if node["type"] == "mol":
            # Check if this is the final product (root node)
            if depth == 0:
                if checker.check_fg("Nitrile", node["smiles"]):
                    nitrile_present_in_final = True
                    print(f"Final product contains nitrile group: {node['smiles']}")
                else:
                    print(f"Final product does NOT contain nitrile group: {node['smiles']}")

        elif node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has nitrile
                product_has_nitrile = checker.check_fg("Nitrile", product)

                # Check if any reactant has nitrile
                reactant_has_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants)

                # If a reactant has nitrile but product doesn't, nitrile wasn't preserved
                if reactant_has_nitrile and not product_has_nitrile:
                    all_steps_preserve_nitrile = False
                    print(f"Nitrile not preserved at depth {depth}, reaction: {rsmi}")

                # If product has nitrile but no reactant does, nitrile was introduced
                if product_has_nitrile and not reactant_has_nitrile:
                    print(f"Nitrile introduced at depth {depth}, reaction: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = nitrile_present_in_final and all_steps_preserve_nitrile
    print(f"Nitrile preservation strategy detected: {result}")
    return result
