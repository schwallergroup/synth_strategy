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
    This function detects a strategy involving multiple ester protections
    with different protecting groups (methyl and tert-butyl).
    """
    methyl_ester_found = False
    tbutyl_ester_found = False

    def dfs_traverse(node, depth=0):
        nonlocal methyl_ester_found, tbutyl_ester_found

        indent = "  " * depth
        print(f"{indent}Traversing node of type: {node['type']}")

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                print(f"{indent}Checking reaction: {rsmi}")

                # Check for esterification reactions (protection)
                if checker.check_fg("Carboxylic acid", reactants) and checker.check_fg(
                    "Ester", product
                ):
                    print(f"{indent}Found potential ester formation from carboxylic acid")

                    # Check for methyl ester
                    if "COC(=O)" in product:
                        print(f"{indent}Found methyl ester")
                        methyl_ester_found = True

                    # Check for tert-butyl ester
                    if "CC(C)(C)OC(=O)" in product or "OC(=O)C(C)(C)C" in product:
                        print(f"{indent}Found tert-butyl ester")
                        tbutyl_ester_found = True

                # Check for ester deprotection reactions
                elif checker.check_fg("Ester", reactants) and checker.check_fg(
                    "Carboxylic acid", product
                ):
                    print(f"{indent}Found potential ester deprotection to carboxylic acid")

                    # Check for methyl ester in reactants
                    if "COC(=O)" in reactants:
                        print(f"{indent}Found methyl ester deprotection")
                        methyl_ester_found = True

                    # Check for tert-butyl ester in reactants
                    if "CC(C)(C)OC(=O)" in reactants or "OC(=O)C(C)(C)C" in reactants:
                        print(f"{indent}Found tert-butyl ester deprotection")
                        tbutyl_ester_found = True

                # Also check for transesterification reactions
                elif checker.check_reaction("Transesterification", rsmi):
                    print(f"{indent}Found transesterification reaction")

                    # Check for methyl ester in product
                    if "COC(=O)" in product:
                        print(f"{indent}Found methyl ester in transesterification")
                        methyl_ester_found = True

                    # Check for tert-butyl ester in product
                    if "CC(C)(C)OC(=O)" in product or "OC(=O)C(C)(C)C" in product:
                        print(f"{indent}Found tert-butyl ester in transesterification")
                        tbutyl_ester_found = True

        # Process molecule nodes to check for esters directly
        elif node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]

            # Check for methyl ester in molecule
            if checker.check_fg("Ester", mol_smiles) and "COC(=O)" in mol_smiles:
                print(f"{indent}Found methyl ester in molecule: {mol_smiles}")
                methyl_ester_found = True

            # Check for tert-butyl ester in molecule
            if checker.check_fg("Ester", mol_smiles) and (
                "CC(C)(C)OC(=O)" in mol_smiles or "OC(=O)C(C)(C)C" in mol_smiles
            ):
                print(f"{indent}Found tert-butyl ester in molecule: {mol_smiles}")
                tbutyl_ester_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Methyl ester found: {methyl_ester_found}")
    print(f"tert-Butyl ester found: {tbutyl_ester_found}")

    # Strategy is present if both types of esters are found
    return methyl_ester_found and tbutyl_ester_found
