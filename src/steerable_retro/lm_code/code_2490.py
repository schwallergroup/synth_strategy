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
    This function detects if the synthetic route involves tetrazole-containing compounds.
    """
    # Flag to track if we found the pattern
    found_pattern = False

    def dfs_traverse(node):
        nonlocal found_pattern

        # Check molecule nodes
        if node["type"] == "mol" and "smiles" in node:
            # Check for tetrazole in molecules
            if checker.check_ring("tetrazole", node["smiles"]):
                print(f"Found tetrazole-containing compound: {node['smiles']}")
                found_pattern = True

        # Check reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                # Check for tetrazole-forming reactions
                rsmi = node["metadata"]["rsmi"]

                # Check specific tetrazole-forming reactions
                tetrazole_reactions = [
                    "tetrazole_terminal",
                    "tetrazole_connect_regioisomere_1",
                    "tetrazole_connect_regioisomere_2",
                    "Azide-nitrile click cycloaddition to tetrazole",
                ]
                if any(checker.check_reaction(rxn, rsmi) for rxn in tetrazole_reactions):
                    print(f"Found specific tetrazole-forming reaction: {rsmi}")
                    found_pattern = True

                # Check if tetrazole is formed in this reaction
                product = rsmi.split(">")[-1]
                if checker.check_ring("tetrazole", product):
                    reactants = rsmi.split(">")[0].split(".")
                    tetrazole_in_reactants = any(
                        checker.check_ring("tetrazole", r) for r in reactants
                    )
                    if not tetrazole_in_reactants:
                        print(f"Found tetrazole-forming reaction: {rsmi}")
                        found_pattern = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_pattern
