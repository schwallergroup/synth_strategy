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
    This function detects if methoxy and fluoro substituents on aromatic rings
    are preserved throughout the synthesis.
    """
    # Track methoxy and fluoro groups by their atom indices in each molecule
    methoxy_groups = {}  # mol_smiles -> list of atom indices
    fluoro_groups = {}  # mol_smiles -> list of atom indices

    # First pass: identify all methoxy and fluoro groups in the route
    def identify_groups(node, depth=0):
        if node["type"] == "mol":
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for methoxy on aromatic rings
                has_methoxy = checker.check_fg("Ether", smiles) and checker.check_ring(
                    "benzene", smiles
                )
                if has_methoxy:
                    # Find specific methoxy groups on aromatic rings
                    methoxy_pattern = Chem.MolFromSmarts("c-O-C")
                    matches = mol.GetSubstructMatches(methoxy_pattern)
                    if matches:
                        methoxy_groups[smiles] = matches

                # Check for fluoro on aromatic rings
                has_fluoro = checker.check_fg("Aromatic halide", smiles) and "F" in smiles
                if has_fluoro:
                    # Find specific fluoro groups on aromatic rings
                    fluoro_pattern = Chem.MolFromSmarts("c-F")
                    matches = mol.GetSubstructMatches(fluoro_pattern)
                    if matches:
                        fluoro_groups[smiles] = matches

                print(f"Molecule at depth {depth}: {smiles}")
                print(f"  Has methoxy: {has_methoxy}, Has fluoro: {has_fluoro}")

        # Process children (reactants in retrosynthesis)
        for child in node.get("children", []):
            identify_groups(child, depth + 1)

    # Second pass: check preservation in each reaction
    preserved_methoxy = True
    preserved_fluoro = True

    def check_preservation(node, depth=0):
        nonlocal preserved_methoxy, preserved_fluoro

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check if product has methoxy/fluoro groups
            product_has_methoxy = checker.check_fg("Ether", product) and checker.check_ring(
                "benzene", product
            )
            product_has_fluoro = checker.check_fg("Aromatic halide", product) and "F" in product

            # Check if any reactant has methoxy/fluoro groups
            reactants_have_methoxy = any(
                checker.check_fg("Ether", r) and checker.check_ring("benzene", r) for r in reactants
            )
            reactants_have_fluoro = any(
                checker.check_fg("Aromatic halide", r) and "F" in r for r in reactants
            )

            # Check preservation
            if reactants_have_methoxy and not product_has_methoxy:
                preserved_methoxy = False
                print(f"  Methoxy not preserved in reaction: {rsmi}")

            if reactants_have_fluoro and not product_has_fluoro:
                preserved_fluoro = False
                print(f"  Fluoro not preserved in reaction: {rsmi}")

        # Process children (reactants in retrosynthesis)
        for child in node.get("children", []):
            check_preservation(child, depth + 1)

    # Check if final product has both groups
    final_product_smiles = route["smiles"]
    final_product = Chem.MolFromSmiles(final_product_smiles)

    final_has_methoxy = checker.check_fg("Ether", final_product_smiles) and checker.check_ring(
        "benzene", final_product_smiles
    )
    final_has_fluoro = (
        checker.check_fg("Aromatic halide", final_product_smiles) and "F" in final_product_smiles
    )

    print(f"Final product: {final_product_smiles}")
    print(f"  Has methoxy: {final_has_methoxy}, Has fluoro: {final_has_fluoro}")

    # Only proceed with preservation check if final product has both groups
    if final_has_methoxy and final_has_fluoro:
        identify_groups(route)
        check_preservation(route)

        result = preserved_methoxy and preserved_fluoro
        print(f"Methoxy and fluoro preservation result: {result}")
        return result
    else:
        print("Final product doesn't have both methoxy and fluoro groups")
        return False
