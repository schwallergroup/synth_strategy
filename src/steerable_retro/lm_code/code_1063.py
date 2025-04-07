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
    This function detects if the synthetic route includes an alcohol to bromide
    transformation as part of the preparation of coupling partners.
    """
    found_transformation = False

    def dfs_traverse(node):
        nonlocal found_transformation

        if node["type"] == "reaction":
            try:
                reaction_smiles = node["metadata"]["rsmi"]
                reactants = reaction_smiles.split(">")[0].split(".")
                product = reaction_smiles.split(">")[-1]

                # Check for alcohol to bromide transformation using specific reaction types
                if checker.check_reaction(
                    "PBr3 and alcohol to alkyl bromide", reaction_smiles
                ) or checker.check_reaction("Appel reaction", reaction_smiles):
                    print(
                        f"Found alcohol to bromide transformation via specific reaction: {reaction_smiles}"
                    )
                    found_transformation = True
                    return

                # Check for alcohol in reactants and bromide in product
                alcohol_in_reactants = False
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                        or checker.check_fg("Aromatic alcohol", reactant)
                    ):
                        alcohol_in_reactants = True
                        print(f"Found alcohol in reactant: {reactant}")
                        break

                if alcohol_in_reactants:
                    # Check for primary, secondary, tertiary halides that contain bromine
                    if (
                        checker.check_fg("Primary halide", product)
                        or checker.check_fg("Secondary halide", product)
                        or checker.check_fg("Tertiary halide", product)
                        or checker.check_fg("Aromatic halide", product)
                    ):

                        # Additional check to ensure it's specifically a bromide
                        mol = Chem.MolFromSmiles(product)
                        if mol:
                            for atom in mol.GetAtoms():
                                if atom.GetSymbol() == "Br":
                                    print(f"Found bromide in product: {product}")
                                    found_transformation = True
                                    return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_transformation
