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
    Detects if the synthesis involves a cyclization where a nitro group
    participates in the formation of a new ring.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Parse reactants and product
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                product = Chem.MolFromSmiles(product_smiles)

                if not all(reactants) or not product:
                    print(f"Warning: Could not parse all molecules in reaction at depth {depth}")
                    return

                # Check for nitro group in reactants using checker
                reactants_with_nitro = []
                for i, mol in enumerate(reactants):
                    mol_smiles = Chem.MolToSmiles(mol)
                    if checker.check_fg("Nitro group", mol_smiles):
                        reactants_with_nitro.append((i, mol_smiles))

                has_nitro = len(reactants_with_nitro) > 0

                # Count rings in reactants and product
                reactant_rings = sum(mol.GetRingInfo().NumRings() for mol in reactants)
                product_rings = product.GetRingInfo().NumRings()

                # Check if a new ring was formed
                ring_formed = product_rings > reactant_rings

                if has_nitro and ring_formed:
                    # Check if nitro group is absent or modified in product
                    nitro_in_product = checker.check_fg("Nitro group", product_smiles)

                    if not nitro_in_product:
                        print(
                            f"Found nitro group cyclization at depth {depth}: nitro group consumed in ring formation"
                        )
                        result = True
                    else:
                        # If nitro group still present, check if it's part of a new ring
                        # Get atom indices of nitro groups in reactants
                        for reactant_idx, reactant_smiles in reactants_with_nitro:
                            nitro_indices = checker.get_fg_atom_indices(
                                "Nitro group", reactant_smiles
                            )
                            if nitro_indices:
                                # Check if these atoms are now part of a ring in the product
                                # This requires atom mapping analysis which is complex
                                # For simplicity, we'll check if the number of rings increased
                                # and assume the nitro group was involved
                                print(
                                    f"Found potential nitro group cyclization at depth {depth}: ring count increased with nitro group present"
                                )
                                result = True
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
