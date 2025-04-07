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
    Detects if the synthetic route builds a complex system with 4+ aromatic/heteroaromatic rings.
    """

    def count_aromatic_rings(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f"Failed to parse SMILES: {smiles}")
            return 0

        # Get ring information
        ring_info = mol.GetRingInfo().AtomRings()
        aromatic_ring_count = 0

        # Common aromatic/heteroaromatic rings to check
        aromatic_rings = [
            "benzene",
            "pyridine",
            "pyrrole",
            "furan",
            "thiophene",
            "imidazole",
            "oxazole",
            "thiazole",
            "pyrazole",
            "isoxazole",
            "isothiazole",
            "triazole",
            "tetrazole",
            "indole",
            "benzofuran",
            "benzothiophene",
            "benzimidazole",
            "benzoxazole",
            "benzothiazole",
            "quinoline",
            "isoquinoline",
            "naphthalene",
            "anthracene",
        ]

        # Count rings that are aromatic
        for ring in ring_info:
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_ring_count += 1
                print(f"Found aromatic ring: {[idx for idx in ring]}")

        # Additional check using the checker functions
        for ring_name in aromatic_rings:
            if checker.check_ring(ring_name, smiles):
                print(f"Found {ring_name} ring in molecule")
                # We don't increment the count here to avoid double counting

        print(f"Total aromatic rings found in {smiles}: {aromatic_ring_count}")
        return aromatic_ring_count

    final_product_rings = 0

    def dfs_traverse(node, depth=0):
        nonlocal final_product_rings

        # Add depth to node for tracking
        node["depth"] = depth

        if node["type"] == "mol" and "smiles" in node and depth == 0:
            # This is the final product
            print(f"Analyzing final product: {node['smiles']}")
            final_product_rings = count_aromatic_rings(node["smiles"])
            print(f"Final product contains {final_product_rings} aromatic rings")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Final count of aromatic rings: {final_product_rings}")
    return final_product_rings >= 4
