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
    Detects a synthesis strategy involving iodination of an aromatic ring.
    """
    found_iodination = False

    def dfs_traverse(node, depth=0):
        nonlocal found_iodination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]

                # Primary check using the reaction checker
                if checker.check_reaction("Aromatic iodination", rsmi):
                    print(f"Found aromatic iodination reaction at depth {depth}: {rsmi}")
                    found_iodination = True
                    return

                # Secondary check: verify iodination manually
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Skip if product doesn't have aromatic halide
                if not checker.check_fg("Aromatic halide", product):
                    pass
                elif "I" in product:  # Only proceed if product contains iodine
                    product_mol = Chem.MolFromSmiles(product)

                    # Check if any reactant contains iodine reagent
                    reactants = reactants_part.split(".")
                    iodine_reagent_present = any("I" in r for r in reactants)

                    if product_mol and iodine_reagent_present:
                        # Check for iodine attached to aromatic ring in product
                        for atom in product_mol.GetAtoms():
                            if atom.GetSymbol() == "I":
                                # Check if iodine is attached to an aromatic atom
                                neighbors = atom.GetNeighbors()
                                if neighbors and any(n.GetIsAromatic() for n in neighbors):
                                    # Verify this is an addition (not already in reactants)
                                    # Check if any reactant already has aromatic iodide
                                    reactant_has_aromatic_iodide = any(
                                        checker.check_fg("Aromatic halide", r) and "I" in r
                                        for r in reactants
                                    )

                                    if not reactant_has_aromatic_iodide:
                                        print(
                                            f"Found aromatic iodination (manual check) at depth {depth}: {rsmi}"
                                        )
                                        found_iodination = True
                                        return
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_iodination
