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
    Detects a synthesis strategy involving a ring opening step.
    """
    has_ring_opening = False

    def dfs_traverse(node, depth=0):
        nonlocal has_ring_opening

        if node["type"] == "reaction":
            # Check if reaction is explicitly marked as ring-breaking
            if node.get("metadata", {}).get("RingBreaker", False):
                print(
                    f"Detected explicitly marked ring opening reaction at depth {depth}"
                )
                has_ring_opening = True
                return

            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")
                product = product_part

                # Convert to RDKit molecules
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r.strip()]

                # Skip if any conversion failed
                if product_mol is None or any(r is None for r in reactant_mols):
                    print(f"Failed to parse SMILES at depth {depth}: {rsmi}")
                    return

                # Count rings in reactants and product
                reactant_ring_count = sum(
                    len(Chem.GetSSSR(mol)) for mol in reactant_mols
                )
                product_ring_count = len(Chem.GetSSSR(product_mol))

                print(
                    f"Ring count at depth {depth}: reactants={reactant_ring_count}, product={product_ring_count}"
                )

                if reactant_ring_count > product_ring_count:
                    print(
                        f"Detected ring opening at depth {depth}: {reactant_ring_count} → {product_ring_count}"
                    )
                    has_ring_opening = True

                    # Additional check for specific ring opening reactions
                    if checker.check_reaction("Retro-Diels-Alder from oxazole", rsmi):
                        print(
                            f"Identified specific ring opening: Retro-Diels-Alder from oxazole"
                        )

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Ring opening strategy detection: {has_ring_opening}")
    return has_ring_opening
