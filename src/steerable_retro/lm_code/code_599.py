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
    This function detects the combined strategy of THP protection and aromatic thioether formation.
    This represents the specific synthetic approach used in this route.
    """
    # Initialize flags for each strategy
    thp_protection_found = False
    thioether_formation_found = False
    hydroxylamine_found = False

    # Helper function to traverse the route and check for strategies
    def traverse_route(node, depth=0):
        nonlocal thp_protection_found, thioether_formation_found, hydroxylamine_found

        # Check molecule nodes for functional groups
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for THP ether (protected alcohol)
            if checker.check_ring("tetrahydropyran", mol_smiles) and checker.check_fg(
                "Ether", mol_smiles
            ):
                print(f"Found THP ether at depth {depth}: {mol_smiles}")
                thp_protection_found = True

            # Check for thioether (aromatic)
            if checker.check_fg("Monosulfide", mol_smiles) and "c" in mol_smiles:
                print(f"Found aromatic thioether at depth {depth}: {mol_smiles}")
                thioether_formation_found = True

            # Check for hydroxylamine derivatives - expanded check
            if "NO" in mol_smiles or "ON" in mol_smiles:
                print(f"Found hydroxylamine derivative at depth {depth}: {mol_smiles}")
                hydroxylamine_found = True

            # Specific check for THP-protected hydroxylamine
            if "NOC1CCCCO1" in mol_smiles or "ONC1CCCCO1" in mol_smiles:
                print(f"Found THP-protected hydroxylamine at depth {depth}: {mol_smiles}")
                hydroxylamine_found = True
                thp_protection_found = True

        # Check reaction nodes for specific transformations
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for THP protection reaction
            if "O1CCCCO1" in rxn_smiles or "OC1CCCCO1" in rxn_smiles:
                print(f"Found THP protection reaction at depth {depth}: {rxn_smiles}")
                thp_protection_found = True

            # Check for thioether formation
            if (
                checker.check_reaction("S-alkylation of thiols", rxn_smiles)
                or checker.check_reaction("S-alkylation of thiols with alcohols", rxn_smiles)
                or (("Sc" in rxn_smiles or "SC" in rxn_smiles) and "c" in rxn_smiles)
            ):
                print(f"Found thioether formation reaction at depth {depth}: {rxn_smiles}")
                thioether_formation_found = True

            # Check for hydroxylamine-related reactions - expanded check
            if "NOH" in rxn_smiles or "NO" in rxn_smiles or "ON" in rxn_smiles:
                print(f"Found hydroxylamine-related reaction at depth {depth}: {rxn_smiles}")
                hydroxylamine_found = True

        # Recursively traverse children
        for child in node.get("children", []):
            traverse_route(child, depth + 1)

    # Start traversal from the root
    traverse_route(route)

    # Check if all three strategies are present
    combined_strategy = thp_protection_found and thioether_formation_found and hydroxylamine_found
    print(f"THP protection found: {thp_protection_found}")
    print(f"Aromatic thioether formation found: {thioether_formation_found}")
    print(f"Hydroxylamine intermediate found: {hydroxylamine_found}")
    print(f"Combined protection-thioether-hydroxylamine strategy present: {combined_strategy}")

    return combined_strategy
