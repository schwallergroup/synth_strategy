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
    This function detects if the synthetic route preserves heterocyclic cores
    (benzofuran and thiophene) throughout the synthesis.
    """
    # Track heterocycles through the synthesis
    preserved_benzofuran = True
    preserved_thiophene = True

    # Check if final product has either heterocycle
    final_product_smiles = route["smiles"]

    # Check for benzofuran
    has_benzofuran_in_product = checker.check_ring("benzofuran", final_product_smiles)
    if not has_benzofuran_in_product:
        # Try to detect benzofuran manually by checking for the pattern
        mol = Chem.MolFromSmiles(final_product_smiles)
        if mol:
            benzofuran_pattern = Chem.MolFromSmarts("c1cccc2c1occ2")
            if benzofuran_pattern:
                has_benzofuran_in_product = (
                    len(mol.GetSubstructMatches(benzofuran_pattern)) > 0
                )

    # Check for thiophene
    has_thiophene_in_product = checker.check_ring("thiophene", final_product_smiles)
    if not has_thiophene_in_product:
        # Check if thiazole is present (contains thiophene-like substructure)
        has_thiophene_in_product = checker.check_ring("thiazole", final_product_smiles)
        if not has_thiophene_in_product:
            mol = Chem.MolFromSmiles(final_product_smiles)
            if mol:
                thiophene_pattern = Chem.MolFromSmarts("c1ccsc1")
                if thiophene_pattern:
                    has_thiophene_in_product = (
                        len(mol.GetSubstructMatches(thiophene_pattern)) > 0
                    )

    print(f"Final product: {final_product_smiles}")
    print(f"Has benzofuran in final product: {has_benzofuran_in_product}")
    print(f"Has thiophene in final product: {has_thiophene_in_product}")

    # If neither heterocycle is in the final product, no need to check preservation
    if not has_benzofuran_in_product and not has_thiophene_in_product:
        return False

    def dfs_traverse(node, depth=0):
        nonlocal preserved_benzofuran, preserved_thiophene

        if node["type"] == "reaction":
            # Get reactants and product from reaction
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if benzofuran is preserved in this reaction
                if has_benzofuran_in_product:
                    product_has_benzofuran = checker.check_ring(
                        "benzofuran", product_smiles
                    )
                    if not product_has_benzofuran:
                        # Try manual detection
                        mol = Chem.MolFromSmiles(product_smiles)
                        if mol:
                            benzofuran_pattern = Chem.MolFromSmarts("c1cccc2c1occ2")
                            if benzofuran_pattern:
                                product_has_benzofuran = (
                                    len(mol.GetSubstructMatches(benzofuran_pattern)) > 0
                                )

                    # Check reactants
                    reactants_have_benzofuran = False
                    for r in reactants_smiles:
                        if checker.check_ring("benzofuran", r):
                            reactants_have_benzofuran = True
                            break
                        # Try manual detection
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            benzofuran_pattern = Chem.MolFromSmarts("c1cccc2c1occ2")
                            if (
                                benzofuran_pattern
                                and len(mol.GetSubstructMatches(benzofuran_pattern)) > 0
                            ):
                                reactants_have_benzofuran = True
                                break

                    # If product has benzofuran but reactants don't, it's not preserved
                    if product_has_benzofuran and not reactants_have_benzofuran:
                        print(f"Benzofuran created in reaction: {rsmi}")
                        preserved_benzofuran = False

                # Check if thiophene is preserved in this reaction
                if has_thiophene_in_product:
                    product_has_thiophene = checker.check_ring(
                        "thiophene", product_smiles
                    )
                    if not product_has_thiophene:
                        # Check for thiazole
                        product_has_thiophene = checker.check_ring(
                            "thiazole", product_smiles
                        )
                        if not product_has_thiophene:
                            # Try manual detection
                            mol = Chem.MolFromSmiles(product_smiles)
                            if mol:
                                thiophene_pattern = Chem.MolFromSmarts("c1ccsc1")
                                if thiophene_pattern:
                                    product_has_thiophene = (
                                        len(mol.GetSubstructMatches(thiophene_pattern))
                                        > 0
                                    )

                    # Check reactants
                    reactants_have_thiophene = False
                    for r in reactants_smiles:
                        if checker.check_ring("thiophene", r):
                            reactants_have_thiophene = True
                            break
                        # Check for thiazole
                        if checker.check_ring("thiazole", r):
                            reactants_have_thiophene = True
                            break
                        # Try manual detection
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            thiophene_pattern = Chem.MolFromSmarts("c1ccsc1")
                            if (
                                thiophene_pattern
                                and len(mol.GetSubstructMatches(thiophene_pattern)) > 0
                            ):
                                reactants_have_thiophene = True
                                break

                    # If product has thiophene but reactants don't, it's not preserved
                    if product_has_thiophene and not reactants_have_thiophene:
                        print(f"Thiophene created in reaction: {rsmi}")
                        preserved_thiophene = False

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return true if at least one heterocycle is in the product and was preserved
    if has_benzofuran_in_product and preserved_benzofuran:
        return True
    if has_thiophene_in_product and preserved_thiophene:
        return True

    return False
