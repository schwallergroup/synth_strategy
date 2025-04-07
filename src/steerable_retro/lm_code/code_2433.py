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
    This function detects if the final step in the synthesis involves a deprotection,
    specifically the removal of a tert-butyl group from an ester.
    """
    final_step_is_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_deprotection

        # The final step is the first reaction we encounter (depth=1)
        if depth == 1 and node["type"] == "reaction":
            print(f"Checking final reaction step at depth {depth}")

            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Reaction SMILES: {rsmi}")
                print(f"Product: {product_smiles}")
                print(f"Reactants: {reactants_smiles}")

                # Check if this is a deprotection reaction
                if (
                    checker.check_reaction("COOH ethyl deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    print("Detected a known deprotection reaction")
                    final_step_is_deprotection = True
                    return

                # Check for tert-butyl ester in reactants and carboxylic acid in product
                for reactant in reactants_smiles:
                    if checker.check_fg("Ester", reactant):
                        print(f"Found ester in reactant: {reactant}")
                        # Check if it's specifically a tert-butyl ester
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            tbutyl_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)")
                            if reactant_mol.HasSubstructMatch(tbutyl_pattern):
                                print("Found tert-butyl ester in reactant")
                                # Check if product has carboxylic acid
                                if checker.check_fg("Carboxylic acid", product_smiles):
                                    print("Found carboxylic acid in product")
                                    final_step_is_deprotection = True
                                    return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root (target molecule)
    dfs_traverse(route)

    return final_step_is_deprotection
