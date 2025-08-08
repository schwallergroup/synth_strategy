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
    Detects a strategy involving the synthesis of an oxygen-rich heterocyclic scaffold
    with multiple oxygen-containing functional groups.
    """
    # Track oxygen-containing functional groups in final product
    oxygen_functional_groups = {
        "ether": 0,
        "ester": 0,
        "alcohol": 0,
        "carboxylic_acid": 0,
        "aldehyde": 0,
        "ketone": 0,
        "phenol": 0,
        "heterocycle": False,
    }

    def dfs_traverse(node, depth=0):
        nonlocal oxygen_functional_groups

        if node["type"] == "mol" and depth == 0:  # Final product
            print(f"Analyzing final product: {node['smiles']}")

            # Check for oxygen-containing heterocycles
            oxygen_heterocycles = [
                "furan",
                "pyran",
                "dioxane",
                "tetrahydrofuran",
                "tetrahydropyran",
                "oxirane",
                "oxetane",
                "oxolane",
                "oxane",
                "dioxolane",
                "dioxolene",
                "trioxane",
                "dioxepane",
            ]

            for ring in oxygen_heterocycles:
                if checker.check_ring(ring, node["smiles"]):
                    oxygen_functional_groups["heterocycle"] = True
                    print(f"Found {ring} in final product")
                    break

            # Check for ethers
            if checker.check_fg("Ether", node["smiles"]):
                ether_matches = checker.get_fg_atom_indices("Ether", node["smiles"])
                oxygen_functional_groups["ether"] = len(ether_matches)
                print(f"Found {oxygen_functional_groups['ether']} ether groups in final product")

            # Check for esters
            if checker.check_fg("Ester", node["smiles"]):
                ester_matches = checker.get_fg_atom_indices("Ester", node["smiles"])
                oxygen_functional_groups["ester"] = len(ester_matches)
                print(f"Found {oxygen_functional_groups['ester']} ester groups in final product")

            # Check for alcohols (various types)
            alcohol_types = [
                "Primary alcohol",
                "Secondary alcohol",
                "Tertiary alcohol",
                "Aromatic alcohol",
                "Enol",
            ]
            for alcohol_type in alcohol_types:
                if checker.check_fg(alcohol_type, node["smiles"]):
                    alcohol_matches = checker.get_fg_atom_indices(alcohol_type, node["smiles"])
                    oxygen_functional_groups["alcohol"] += len(alcohol_matches)
            print(f"Found {oxygen_functional_groups['alcohol']} alcohol groups in final product")

            # Check for carboxylic acids
            if checker.check_fg("Carboxylic acid", node["smiles"]):
                acid_matches = checker.get_fg_atom_indices("Carboxylic acid", node["smiles"])
                oxygen_functional_groups["carboxylic_acid"] = len(acid_matches)
                print(
                    f"Found {oxygen_functional_groups['carboxylic_acid']} carboxylic acid groups in final product"
                )

            # Check for aldehydes
            if checker.check_fg("Aldehyde", node["smiles"]):
                aldehyde_matches = checker.get_fg_atom_indices("Aldehyde", node["smiles"])
                oxygen_functional_groups["aldehyde"] = len(aldehyde_matches)
                print(
                    f"Found {oxygen_functional_groups['aldehyde']} aldehyde groups in final product"
                )

            # Check for ketones
            if checker.check_fg("Ketone", node["smiles"]):
                ketone_matches = checker.get_fg_atom_indices("Ketone", node["smiles"])
                oxygen_functional_groups["ketone"] = len(ketone_matches)
                print(f"Found {oxygen_functional_groups['ketone']} ketone groups in final product")

            # Check for phenols
            if checker.check_fg("Phenol", node["smiles"]):
                phenol_matches = checker.get_fg_atom_indices("Phenol", node["smiles"])
                oxygen_functional_groups["phenol"] = len(phenol_matches)
                print(f"Found {oxygen_functional_groups['phenol']} phenol groups in final product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if the strategy is present
    total_o_groups = (
        oxygen_functional_groups["ether"]
        + oxygen_functional_groups["ester"]
        + oxygen_functional_groups["alcohol"]
        + oxygen_functional_groups["carboxylic_acid"]
        + oxygen_functional_groups["aldehyde"]
        + oxygen_functional_groups["ketone"]
        + oxygen_functional_groups["phenol"]
    )

    strategy_present = (
        oxygen_functional_groups["heterocycle"]
        and total_o_groups
        >= 2  # At least 2 additional oxygen functional groups besides the heterocycle
    )

    print(f"Strategy detection result: {strategy_present}")
    print(f"- Oxygen heterocycle: {oxygen_functional_groups['heterocycle']}")
    print(f"- Total oxygen functional groups: {total_o_groups}")
    print(f"- Functional group breakdown: {oxygen_functional_groups}")

    return strategy_present
