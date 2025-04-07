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
    Detects if the synthesis route involves ester hydrolysis to carboxylic acid.
    """
    has_ester_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal has_ester_hydrolysis

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Depth {depth}: Analyzing reaction: {rsmi}")

                # Check if this is an ester hydrolysis reaction
                if checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                ):
                    print(f"Depth {depth}: Found ester hydrolysis reaction")

                    # Verify ester in reactants and carboxylic acid in product
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol is not None:
                        has_ester = any(
                            r is not None and checker.check_fg("Ester", Chem.MolToSmiles(r))
                            for r in reactant_mols
                        )
                        has_acid = checker.check_fg("Carboxylic acid", product_smiles)

                        if has_ester and has_acid:
                            print(f"Depth {depth}: Confirmed ester hydrolysis to carboxylic acid")
                            has_ester_hydrolysis = True

                # Alternative check using functional groups if reaction check fails
                if not has_ester_hydrolysis:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol is not None:
                        has_ester = any(
                            r is not None and checker.check_fg("Ester", Chem.MolToSmiles(r))
                            for r in reactant_mols
                        )
                        has_acid = checker.check_fg("Carboxylic acid", product_smiles)

                        if (
                            has_ester
                            and has_acid
                            and not any(
                                checker.check_fg("Carboxylic acid", Chem.MolToSmiles(r))
                                for r in reactant_mols
                                if r is not None
                            )
                        ):
                            print(
                                f"Depth {depth}: Detected ester hydrolysis based on functional groups"
                            )
                            has_ester_hydrolysis = True

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    print("Starting analysis of synthesis route for ester hydrolysis")
    dfs_traverse(route)
    print(f"Ester hydrolysis detected: {has_ester_hydrolysis}")
    return has_ester_hydrolysis
