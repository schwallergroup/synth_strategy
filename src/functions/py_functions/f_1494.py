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
    This function detects if an aromatic fluorine is retained throughout the synthesis.
    """
    # Track if we've found at least one aromatic fluorine that's retained
    aromatic_fluorine_retained = False

    def has_aromatic_fluorine(smiles):
        """Check if a molecule has an aromatic carbon with fluorine attached"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Check for aromatic halide
        if checker.check_fg("Aromatic halide", smiles):
            # Verify it's specifically fluorine
            for atom in mol.GetAtoms():
                # Check if it's fluorine
                if atom.GetAtomicNum() == 9:
                    # Check if it's connected to an aromatic carbon
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIsAromatic():
                            return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_fluorine_retained

        if node["type"] == "mol":
            # Check if this molecule has aromatic fluorine
            if has_aromatic_fluorine(node["smiles"]):
                print(
                    f"Found molecule with aromatic fluorine at depth {depth}: {node['smiles']}"
                )

                # If this is the target molecule (depth 0), mark as retained
                if depth == 0:
                    aromatic_fluorine_retained = True

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has aromatic fluorine
            product_has_aromatic_f = has_aromatic_fluorine(product)

            # Check if any reactant has aromatic fluorine
            reactants_have_aromatic_f = any(has_aromatic_fluorine(r) for r in reactants)

            # If product has aromatic fluorine but reactants don't, it was introduced in this step
            if product_has_aromatic_f and not reactants_have_aromatic_f:
                print(f"Aromatic fluorine introduced at reaction: {rsmi}")

            # If reactants have aromatic fluorine but product doesn't, it was removed in this step
            if reactants_have_aromatic_f and not product_has_aromatic_f:
                print(f"Aromatic fluorine removed at reaction: {rsmi}")
                aromatic_fluorine_retained = False

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if aromatic_fluorine_retained:
        print("Aromatic fluorine is retained throughout the synthesis")
    else:
        print("Aromatic fluorine is NOT retained throughout the synthesis")

    return aromatic_fluorine_retained
