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
    Detects a strategy involving the synthesis of a compound containing multiple
    nitrogen heterocycles (pyrazole, thiazole, isoxazole, tetrazole, etc.).
    """
    # Track nitrogen heterocycles in the final product
    nitrogen_heterocycles_in_final = set()

    # Track which heterocycles are formed during synthesis
    heterocycle_formation = {
        "pyrazole": False,
        "thiazole": False,
        "isoxazole": False,
        "tetrazole": False,
        "triazole": False,
        "imidazole": False,
    }

    def dfs_traverse(node, depth=0):
        nonlocal nitrogen_heterocycles_in_final

        # Check if this is a molecule node
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if this is the root node (final product)
            if depth == 0:
                # Check for nitrogen heterocycles in the final product
                for ring_name in [
                    "pyrazole",
                    "thiazole",
                    "isoxazole",
                    "tetrazole",
                    "triazole",
                    "imidazole",
                ]:
                    if checker.check_ring(ring_name, mol_smiles):
                        nitrogen_heterocycles_in_final.add(ring_name)
                        print(f"Final product contains {ring_name}: {mol_smiles}")

                # Manual check for isoxazole pattern if the checker didn't detect it
                if "isoxazole" not in nitrogen_heterocycles_in_final and "=NOC" in mol_smiles:
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        for atom in mol.GetAtoms():
                            if atom.GetSymbol() == "N" and atom.GetIsAromatic():
                                for neighbor in atom.GetNeighbors():
                                    if neighbor.GetSymbol() == "O" and neighbor.GetIsAromatic():
                                        nitrogen_heterocycles_in_final.add("isoxazole")
                                        print(
                                            f"Manually detected isoxazole in final product: {mol_smiles}"
                                        )
                                        break

                print(f"Final product SMILES: {mol_smiles}")
                print(f"Nitrogen heterocycles in final product: {nitrogen_heterocycles_in_final}")

        # Check if this is a reaction node
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for specific heterocycle formation reactions
            if checker.check_reaction("pyrazole", rsmi):
                print(f"Pyrazole formation reaction detected: {rsmi}")
                heterocycle_formation["pyrazole"] = True

            if checker.check_reaction("thiazole", rsmi):
                print(f"Thiazole formation reaction detected: {rsmi}")
                heterocycle_formation["thiazole"] = True

            if checker.check_reaction("{oxadiazole}", rsmi) or checker.check_reaction(
                "isoxazole", rsmi
            ):
                print(f"Isoxazole/oxadiazole formation reaction detected: {rsmi}")
                heterocycle_formation["isoxazole"] = True

            if (
                checker.check_reaction("{tetrazole_terminal}", rsmi)
                or checker.check_reaction("{tetrazole_connect_regioisomere_1}", rsmi)
                or checker.check_reaction("{tetrazole_connect_regioisomere_2}", rsmi)
            ):
                print(f"Tetrazole formation reaction detected: {rsmi}")
                heterocycle_formation["tetrazole"] = True

            if checker.check_reaction(
                "{1,2,4-triazole_acetohydrazide}", rsmi
            ) or checker.check_reaction("{1,2,4-triazole_carboxylic-acid/ester}", rsmi):
                print(f"Triazole formation reaction detected: {rsmi}")
                heterocycle_formation["triazole"] = True

            # Check if heterocycles appear in product but not in reactants
            for ring_name in heterocycle_formation.keys():
                any_reactant_has_ring = any(checker.check_ring(ring_name, r) for r in reactants)
                product_has_ring = checker.check_ring(ring_name, product)

                if not any_reactant_has_ring and product_has_ring:
                    print(f"{ring_name.capitalize()} formation detected in reaction: {rsmi}")
                    heterocycle_formation[ring_name] = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Print final status for debugging
    print(f"Nitrogen heterocycles in final product: {nitrogen_heterocycles_in_final}")
    print(f"Heterocycle formation during synthesis: {heterocycle_formation}")

    # Return True if the final product contains at least two nitrogen heterocycles
    # and at least one of them is formed during synthesis
    has_multiple_n_heterocycles = len(nitrogen_heterocycles_in_final) >= 2
    any_heterocycle_formed = any(heterocycle_formation.values())

    return has_multiple_n_heterocycles and any_heterocycle_formed
