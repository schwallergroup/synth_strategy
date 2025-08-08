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
    Detects if the route contains a Williamson ether synthesis
    (alcohol + alkyl halide → ether).
    """
    ether_synthesis_found = False

    def dfs_traverse(node):
        nonlocal ether_synthesis_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]

            # First, directly check if this is a Williamson ether synthesis reaction
            if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                print(f"Found Williamson ether synthesis via reaction check: {rsmi}")
                ether_synthesis_found = True
                return

            # If direct check fails, perform more detailed analysis
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol in reactants
            alcohol_present = any(
                checker.check_fg("Primary alcohol", r)
                or checker.check_fg("Secondary alcohol", r)
                or checker.check_fg("Tertiary alcohol", r)
                or checker.check_fg("Phenol", r)
                for r in reactants
            )

            # Check for alkyl halide in reactants
            halide_present = any(
                checker.check_fg("Primary halide", r)
                or checker.check_fg("Secondary halide", r)
                or checker.check_fg("Tertiary halide", r)
                for r in reactants
            )

            # Check for ether in product
            ether_formed = checker.check_fg("Ether", product)

            if alcohol_present and halide_present and ether_formed:
                # Additional verification that this is indeed a Williamson ether synthesis
                # by checking the reaction pattern more broadly
                if "O" in product and any("Br" in r or "Cl" in r or "I" in r for r in reactants):
                    print(f"Found Williamson ether synthesis via functional group analysis: {rsmi}")
                    ether_synthesis_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ether_synthesis_found
