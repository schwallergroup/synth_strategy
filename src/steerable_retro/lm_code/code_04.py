#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects a linear synthesis strategy involving sequential functionalization
    of a heterocyclic core.
    """
    # List of common heterocycles to check
    heterocycle_rings = [
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "purine",
    ]

    # Track the synthesis path and heterocycle modifications
    linear_path = []
    heterocycle_modifications = []

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        current_path = path + [node]

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if this molecule contains a heterocycle
            heterocycle_found = False
            for ring_name in heterocycle_rings:
                if checker.check_ring(ring_name, mol_smiles):
                    heterocycle_found = True
                    print(f"Found heterocycle {ring_name} in molecule: {mol_smiles}")
                    break

            # If this is a leaf node (starting material) with a heterocycle, record it
            if heterocycle_found and node.get("in_stock", False):
                heterocycle_modifications.append((depth, mol_smiles, "starting_material"))

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains a heterocycle
                product_heterocycle = None
                for ring_name in heterocycle_rings:
                    if checker.check_ring(ring_name, product):
                        product_heterocycle = ring_name
                        break

                # Check if any reactant contains a heterocycle
                reactant_heterocycle = None
                for reactant in reactants:
                    for ring_name in heterocycle_rings:
                        if checker.check_ring(ring_name, reactant):
                            reactant_heterocycle = ring_name
                            break
                    if reactant_heterocycle:
                        break

                # If both reactant and product have heterocycles, this is a modification step
                if reactant_heterocycle and product_heterocycle:
                    print(
                        f"Found heterocycle modification at depth {depth}: {reactant_heterocycle} -> {product_heterocycle}"
                    )
                    heterocycle_modifications.append((depth, rsmi, "modification"))
                    linear_path.append((depth, current_path))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # Sort modifications by depth (retrosynthetic order)
    heterocycle_modifications.sort(key=lambda x: x[0])
    linear_path.sort(key=lambda x: x[0])

    print(f"Heterocycle modifications: {heterocycle_modifications}")
    print(f"Linear path: {[depth for depth, _ in linear_path]}")

    # Check if we have heterocycle modifications
    if len(heterocycle_modifications) < 2:
        print("Not enough heterocycle modifications found")
        return False

    # Check if the modifications form a linear sequence (increasing depths, not necessarily consecutive)
    is_linear = True
    if len(linear_path) > 1:
        depths = [depth for depth, _ in linear_path]
        for i in range(1, len(depths)):
            if depths[i] <= depths[i - 1]:  # Check if depths are strictly increasing
                is_linear = False
                print(f"Non-increasing depths: {depths[i-1]} and {depths[i]}")
                break
    else:
        is_linear = False
        print("Not enough reactions in path")

    # Check if the same heterocycle core is being modified throughout
    same_core = True
    if len(heterocycle_modifications) > 1:
        # Extract the heterocycle types from each modification
        heterocycle_types = []
        for _, smiles, mod_type in heterocycle_modifications:
            if mod_type == "starting_material":
                for ring_name in heterocycle_rings:
                    if checker.check_ring(ring_name, smiles):
                        heterocycle_types.append(ring_name)
                        print(f"Starting material heterocycle: {ring_name}")
                        break
            else:  # modification
                product = smiles.split(">")[-1]
                for ring_name in heterocycle_rings:
                    if checker.check_ring(ring_name, product):
                        heterocycle_types.append(ring_name)
                        print(f"Modification product heterocycle: {ring_name}")
                        break

        # Check if all modifications involve the same heterocycle type
        if len(set(heterocycle_types)) > 1:
            same_core = False
            print(f"Multiple heterocycle types found: {set(heterocycle_types)}")

    print(f"Is linear: {is_linear}")
    print(f"Same heterocycle core: {same_core}")

    # Return True if it's a linear synthesis with consistent heterocycle modifications
    return is_linear and same_core and len(heterocycle_modifications) >= 2
