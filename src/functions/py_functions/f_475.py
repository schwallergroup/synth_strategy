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
    Detects if the synthesis route includes an ester reduction to alcohol step.
    """
    found_ester_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ester_reduction

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            try:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")
                product = product_part

                # Check if this is a reduction of ester to primary alcohol reaction
                if checker.check_reaction(
                    "Reduction of ester to primary alcohol", rsmi
                ):
                    print(f"Found ester reduction to alcohol at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    found_ester_reduction = True
                    return

                # Fallback check: verify ester in reactants and primary alcohol in product
                reactant_has_ester = any(
                    checker.check_fg("Ester", reactant) for reactant in reactants
                )
                product_has_alcohol = checker.check_fg("Primary alcohol", product)

                if reactant_has_ester and product_has_alcohol:
                    # Additional check to ensure it's a reduction reaction
                    # Look for patterns consistent with ester reduction
                    reactant_mol = Chem.MolFromSmiles(reactants[0])
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        # Check if carbon count is preserved (no carbon-carbon bond breaking)
                        if reactant_mol.GetNumAtoms(
                            onlyExplicit=True
                        ) >= product_mol.GetNumAtoms(onlyExplicit=True):
                            print(
                                f"Found potential ester reduction to alcohol at depth {depth}"
                            )
                            print(f"Reaction SMILES: {rsmi}")
                            found_ester_reduction = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_ester_reduction
