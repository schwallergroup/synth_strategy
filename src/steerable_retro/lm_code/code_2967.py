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
    Detects if the synthesis involves protection of a terminal alkyne with a TMS group.
    """
    tms_protection = False

    def dfs_traverse(node):
        nonlocal tms_protection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for terminal alkyne in reactants and TMS group in product
            terminal_alkyne_reactant = False
            for reactant in reactants:
                # Check if reactant has a terminal alkyne
                if checker.check_fg("Alkyne", reactant):
                    # Verify it's a terminal alkyne by checking the molecule
                    try:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            terminal_alkyne_pattern = Chem.MolFromSmarts("C#[CH]")
                            if r_mol.HasSubstructMatch(terminal_alkyne_pattern):
                                terminal_alkyne_reactant = True
                                print(f"Terminal alkyne found in reactant: {reactant}")
                    except Exception as e:
                        print(f"Error checking terminal alkyne: {e}")

            # Check for TMS-protected alkyne in product
            tms_alkyne_product = False
            try:
                p_mol = Chem.MolFromSmiles(product)
                if p_mol:
                    tms_pattern = Chem.MolFromSmarts("C#C[Si](C)(C)C")
                    if p_mol.HasSubstructMatch(tms_pattern):
                        tms_alkyne_product = True
                        print(f"TMS-protected alkyne found in product: {product}")
            except Exception as e:
                print(f"Error checking TMS-alkyne: {e}")

            # If we found both a terminal alkyne in reactants and a TMS-alkyne in product,
            # this is likely a TMS protection reaction
            if terminal_alkyne_reactant and tms_alkyne_product:
                tms_protection = True
                print(f"TMS-alkyne protection detected in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return tms_protection
