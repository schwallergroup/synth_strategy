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
    This function detects synthesis routes involving trifluoromethyl-containing compounds.
    It checks both molecule nodes and reaction nodes for the presence of trifluoromethyl groups.
    """
    has_trifluoromethyl = False

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoromethyl

        # Process molecule nodes
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]
                if checker.check_fg("Trifluoro group", mol_smiles):
                    print(f"Detected trifluoromethyl group in molecule: {mol_smiles}")
                    has_trifluoromethyl = True
            except Exception as e:
                print(f"Error processing molecule node: {e}")

        # Process reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant contains trifluoromethyl group
                for reactant in reactants:
                    if checker.check_fg("Trifluoro group", reactant):
                        print(f"Detected trifluoromethyl group in reactant: {reactant}")
                        has_trifluoromethyl = True

                # Check if product contains trifluoromethyl group
                if checker.check_fg("Trifluoro group", product):
                    print(f"Detected trifluoromethyl group in product: {product}")
                    has_trifluoromethyl = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_trifluoromethyl
