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
    Detects a synthesis strategy involving the preparation of a thioether-containing indole derivative
    that is later used in a coupling reaction.
    """
    # Initialize tracking variables
    has_thioether = False
    has_indole = False
    has_thioether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_thioether, has_indole, has_thioether_formation

        if node["type"] == "mol":
            if "smiles" in node:
                # Check for indole
                if checker.check_ring("indole", node["smiles"]):
                    has_indole = True
                    print(f"Detected indole ring in molecule: {node['smiles']}")

                # Check for thioether (monosulfide)
                if checker.check_fg("Monosulfide", node["smiles"]):
                    has_thioether = True
                    print(f"Detected thioether group in molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for thioether formation reactions
            if (
                checker.check_reaction("S-alkylation of thiols", rsmi)
                or checker.check_reaction("S-alkylation of thiols (ethyl)", rsmi)
                or checker.check_reaction("thioether_nucl_sub", rsmi)
            ):
                has_thioether_formation = True
                print(f"Detected thioether formation reaction: {rsmi}")

            # Additional check for other potential thioether formation reactions
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has a thiol group and product has a thioether
            thiol_in_reactants = any(
                checker.check_fg("Aliphatic thiol", r) or checker.check_fg("Aromatic thiol", r)
                for r in reactants
            )

            if (
                thiol_in_reactants
                and checker.check_fg("Monosulfide", product)
                and not has_thioether_formation
            ):
                has_thioether_formation = True
                print(f"Detected thioether formation from thiol: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all required elements are present
    result = has_thioether and has_indole and has_thioether_formation

    print(f"Thioether detected: {has_thioether}")
    print(f"Indole detected: {has_indole}")
    print(f"Thioether formation detected: {has_thioether_formation}")
    print(f"Overall thioether-indole strategy detected: {result}")

    return result
