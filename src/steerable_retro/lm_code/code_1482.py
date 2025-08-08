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
    Detects if the synthesis route preserves a nitrile group from starting materials to final product.
    """
    # First, check if the final product contains a nitrile group
    if route["type"] != "mol" or not checker.check_fg("Nitrile", route["smiles"]):
        return False

    # Track if we find a nitrile in any starting material
    nitrile_in_starting_material = False

    def dfs(node, depth=0):
        nonlocal nitrile_in_starting_material

        if node["type"] == "mol" and node.get("in_stock", False) and node["smiles"]:
            if checker.check_fg("Nitrile", node["smiles"]):
                nitrile_in_starting_material = True
                print(f"Found nitrile in starting material: {node['smiles']}")

        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)

    # Check if there are no reactions that would typically modify a nitrile
    nitrile_modified = False

    def check_nitrile_modification(node, depth=0):
        nonlocal nitrile_modified

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for reactions that typically modify nitriles
            if (
                checker.check_reaction("Nitrile to amide", rxn_smiles)
                or checker.check_reaction("Reduction of nitrile to amine", rxn_smiles)
                or checker.check_reaction("Oxidation of nitrile to carboxylic acid", rxn_smiles)
            ):
                nitrile_modified = True
                print(f"Found nitrile-modifying reaction: {rxn_smiles}")
            else:
                # Additional check: see if nitrile is lost in any reaction
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                nitrile_in_reactants = any(checker.check_fg("Nitrile", r) for r in reactants)
                nitrile_in_product = checker.check_fg("Nitrile", product)

                if nitrile_in_reactants and not nitrile_in_product:
                    nitrile_modified = True
                    print(f"Found nitrile-modifying reaction (by FG analysis): {rxn_smiles}")

        for child in node.get("children", []):
            check_nitrile_modification(child, depth + 1)

    check_nitrile_modification(route)

    return nitrile_in_starting_material and not nitrile_modified
