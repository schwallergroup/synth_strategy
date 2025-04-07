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
    Detects late-stage heterocycle formation involving a nitro group.
    Looks for cyclization reactions in the final steps (depth 0-1) where
    a nitro group is involved in the formation of a new ring.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and depth <= 1:
            # Check if this is a late-stage reaction (depth 0 or 1)
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Parse reactants and product
                reactant_list = reactants_smiles.split(".")
                reactants = [Chem.MolFromSmiles(r) for r in reactant_list if r]
                product = Chem.MolFromSmiles(product_smiles)

                if not all(reactants) or not product:
                    print(f"Warning: Could not parse all molecules in reaction at depth {depth}")
                    return

                # Check for nitro group in reactants using checker function
                has_nitro = False
                nitro_reactant_idx = -1
                for i, r_smiles in enumerate(reactant_list):
                    if checker.check_fg("Nitro group", r_smiles):
                        has_nitro = True
                        nitro_reactant_idx = i
                        print(f"Found nitro group in reactant at depth {depth}")
                        break

                if has_nitro:
                    # Count rings in reactants and product
                    reactant_rings = sum(mol.GetRingInfo().NumRings() for mol in reactants)
                    product_rings = product.GetRingInfo().NumRings()

                    # Check if a new ring was formed
                    if product_rings > reactant_rings:
                        print(f"Ring count increased: {reactant_rings} -> {product_rings}")

                        # Check if the product contains a heterocycle
                        heterocycle_found = False
                        for ring_name in [
                            "pyrrole",
                            "pyridine",
                            "pyrazole",
                            "imidazole",
                            "oxazole",
                            "thiazole",
                            "pyrimidine",
                            "pyrazine",
                            "triazole",
                            "tetrazole",
                            "indole",
                            "quinoline",
                            "isoquinoline",
                            "benzimidazole",
                            "benzoxazole",
                            "benzothiazole",
                        ]:
                            if checker.check_ring(ring_name, product_smiles):
                                heterocycle_found = True
                                print(f"Found heterocycle: {ring_name} in product at depth {depth}")
                                break

                        if heterocycle_found:
                            # Check if the nitro group is involved in the ring formation
                            # This is a simplification - in reality, the nitro group might be
                            # reduced to an amine first, which then participates in ring formation
                            if not checker.check_fg("Nitro group", product_smiles):
                                print(
                                    f"Nitro group was transformed during heterocycle formation at depth {depth}"
                                )
                                result = True
                            elif checker.check_fg("Nitro group", product_smiles):
                                # If nitro group is still present, check if it's part of the heterocycle
                                # This is a simplified check - ideally we would check atom mappings
                                print(
                                    f"Nitro group still present in product, may be part of heterocycle at depth {depth}"
                                )
                                result = True
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
