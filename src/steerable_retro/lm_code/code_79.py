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
    Detects a synthetic strategy where modifications are made to substituents on a pyridine ring,
    while maintaining the core pyridine scaffold throughout the synthesis.
    """
    # Initialize tracking variables
    main_intermediates = []

    def identify_main_pathway(node, depth=0):
        """Identify the main synthetic pathway molecules"""
        if node["type"] == "mol" and not node.get("in_stock", False):
            # Add target molecule (depth 0) to main intermediates
            if depth == 0:
                main_intermediates.append(node["smiles"])
                print(f"Added target molecule to main pathway: {node['smiles']}")

            # For reaction nodes, identify the main product
            for child in node.get("children", []):
                if (
                    child["type"] == "reaction"
                    and "metadata" in child
                    and "rsmi" in child["metadata"]
                ):
                    rsmi = child["metadata"]["rsmi"]
                    product = rsmi.split(">")[-1]

                    # The product of this reaction should be the current molecule
                    if product == node["smiles"]:
                        # Find the main reactant (usually the first/largest one)
                        reactants = rsmi.split(">")[0].split(".")

                        # Filter out small molecules (likely reagents)
                        main_reactants = [r for r in reactants if len(r) > 15]

                        if main_reactants:
                            # Add all significant reactants to main pathway
                            for reactant in main_reactants:
                                main_intermediates.append(reactant)
                                print(f"Added reactant to main pathway: {reactant}")

        # Traverse children
        for child in node.get("children", []):
            identify_main_pathway(child, depth + 1)

    # Start traversal to identify main pathway
    identify_main_pathway(route)

    # Check if all main intermediates contain pyridine
    pyridine_count = 0
    for mol_smiles in main_intermediates:
        if checker.check_ring("pyridine", mol_smiles):
            pyridine_count += 1
            print(f"Main intermediate with pyridine: {mol_smiles}")
        else:
            print(f"Main intermediate WITHOUT pyridine: {mol_smiles}")

    # Ensure we have at least one pyridine-containing molecule and all have pyridine
    has_pyridine_strategy = pyridine_count > 0 and pyridine_count == len(main_intermediates)

    print(f"Pyridine molecules in main pathway: {pyridine_count}/{len(main_intermediates)}")
    print(f"Pyridine modification strategy detected: {has_pyridine_strategy}")

    return has_pyridine_strategy
