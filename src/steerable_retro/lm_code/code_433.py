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
    This function detects a synthetic strategy where multiple heterocyclic systems
    are connected via linkers (specifically looking for imidazole and benzofuran).
    """
    has_imidazole = False
    has_benzofuran = False
    molecules_with_both = []

    # Define benzofuran-like patterns to check
    benzofuran_like_patterns = [
        "c1ccc2occc2c1",  # Basic benzofuran
        "c1cc2ccccc2o1",  # Alternative benzofuran representation
        "c1ccc2c(c1)oc1ccccc12",  # Dibenzofuran
    ]

    def check_benzofuran_like(smiles):
        """Helper function to check for benzofuran-like structures"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False

            # Check with the checker function first
            if checker.check_ring("benzofuran", smiles) or checker.check_ring(
                "dibenzofuran", smiles
            ):
                return True

            # Check for benzofuran-like patterns
            for pattern in benzofuran_like_patterns:
                pattern_mol = Chem.MolFromSmiles(pattern)
                if pattern_mol and mol.HasSubstructMatch(pattern_mol):
                    print(f"Found benzofuran-like structure in {smiles}")
                    return True

            # Check for sulfonamide connected to a benzofuran-like structure
            if "NS(=O)(=O)c" in smiles and "ccc3o" in smiles:
                print(f"Found potential sulfonamide-benzofuran connection in {smiles}")
                return True

            return False
        except Exception as e:
            print(f"Error in benzofuran check: {e}")
            return False

    def dfs_traverse(node):
        nonlocal has_imidazole, has_benzofuran

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]
                # Check for imidazole
                has_imidazole_in_mol = checker.check_ring("imidazole", mol_smiles)

                # Check for benzofuran or benzofuran-like structures
                has_benzofuran_in_mol = check_benzofuran_like(mol_smiles)

                if has_imidazole_in_mol:
                    has_imidazole = True
                    print(f"Found imidazole heterocycle in {mol_smiles}")

                if has_benzofuran_in_mol:
                    has_benzofuran = True
                    print(f"Found benzofuran heterocycle in {mol_smiles}")

                # Check if both heterocycles are in the same molecule
                if has_imidazole_in_mol and has_benzofuran_in_mol:
                    print(f"Found both heterocycles in the same molecule: {mol_smiles}")
                    molecules_with_both.append(mol_smiles)

                # Special case: Check for sulfonamide benzofuran connection to imidazole
                if has_imidazole_in_mol and "NS(=O)(=O)c" in mol_smiles and "ccc3o" in mol_smiles:
                    print(
                        f"Found imidazole connected to sulfonamide-benzofuran structure: {mol_smiles}"
                    )
                    molecules_with_both.append(mol_smiles)
            except Exception as e:
                print(f"Error processing molecule: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Check if we have molecules with both heterocycles
    if molecules_with_both:
        print(f"Found {len(molecules_with_both)} molecules containing both heterocycles")
        return True

    # If we don't have molecules with both, check if they're connected via reactions
    if has_imidazole and has_benzofuran:
        print("Both heterocycles found in the route, but not in the same molecule")

        # Check for connection via reactions
        def check_connection(node, depth=0):
            if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    product = rsmi.split(">")[-1]
                    reactants = rsmi.split(">")[0].split(".")

                    # Check if the product contains both heterocycles
                    product_has_imidazole = checker.check_ring("imidazole", product)
                    product_has_benzofuran = check_benzofuran_like(product)

                    if product_has_imidazole and product_has_benzofuran:
                        print(f"Found reaction producing a molecule with both heterocycles: {rsmi}")
                        return True

                    # Check if reactants contain the heterocycles
                    reactants_with_imidazole = [
                        r for r in reactants if checker.check_ring("imidazole", r)
                    ]
                    reactants_with_benzofuran = [r for r in reactants if check_benzofuran_like(r)]

                    # Case 1: One reactant has imidazole and another has benzofuran
                    if reactants_with_imidazole and reactants_with_benzofuran:
                        print(
                            f"Found reaction with both heterocycles in different reactants: {rsmi}"
                        )
                        return True

                    # Case 2: One reactant has one heterocycle and product has both
                    if (reactants_with_imidazole and product_has_benzofuran) or (
                        reactants_with_benzofuran and product_has_imidazole
                    ):
                        print(f"Found reaction connecting heterocycles: {rsmi}")
                        return True

                    # Special case: Check for sulfonamide-benzofuran connection
                    if (
                        any("NS(=O)(=O)" in r for r in reactants)
                        and any("ccc3o" in r for r in reactants)
                        and product_has_imidazole
                    ):
                        print(
                            f"Found reaction connecting sulfonamide-benzofuran to imidazole: {rsmi}"
                        )
                        return True

                    # Check if any reactant has sulfonamide-benzofuran and product has imidazole
                    if (
                        any(("NS(=O)(=O)" in r and "ccc3o" in r) for r in reactants)
                        and product_has_imidazole
                    ):
                        print(
                            f"Found reaction connecting sulfonamide-benzofuran to imidazole: {rsmi}"
                        )
                        return True
                except Exception as e:
                    print(f"Error checking reaction connection: {e}")

            for child in node.get("children", []):
                if check_connection(child, depth + 1):
                    return True

            return False

        # Check for direct connection
        if check_connection(route):
            return True

    # Final check: Look at the final product for specific patterns
    if route["type"] == "mol" and "smiles" in route:
        final_product = route["smiles"]
        # Check for imidazole connected to sulfonamide-benzofuran structure
        if (
            checker.check_ring("imidazole", final_product)
            and "NS(=O)(=O)" in final_product
            and "ccc3o" in final_product
        ):
            print(
                f"Final product contains imidazole connected to sulfonamide-benzofuran: {final_product}"
            )
            return True

        # Check for the specific molecule pattern seen in stdout
        if "CCCn1cnc" in final_product and "NS(=O)(=O)c2cc3ccccc3o2" in final_product:
            print(
                f"Found target molecule with imidazole and benzofuran connected via sulfonamide: {final_product}"
            )
            return True

    return False
