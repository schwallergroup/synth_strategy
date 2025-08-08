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
    Detects if there's a strategic bromination step in the route.
    """
    bromination_found = False

    def dfs(node):
        nonlocal bromination_found

        # Check if this is a reaction node
        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a bromination reaction
            is_bromination = (
                checker.check_reaction("Aromatic bromination", rsmi)
                or checker.check_reaction("Bromination", rsmi)
                or checker.check_reaction("Wohl-Ziegler bromination benzyl primary", rsmi)
                or checker.check_reaction("Wohl-Ziegler bromination benzyl secondary", rsmi)
                or checker.check_reaction("Wohl-Ziegler bromination benzyl tertiary", rsmi)
                or checker.check_reaction("Wohl-Ziegler bromination allyl primary", rsmi)
                or checker.check_reaction("Wohl-Ziegler bromination allyl secondary", rsmi)
                or checker.check_reaction("Wohl-Ziegler bromination allyl tertiary", rsmi)
            )

            # Also check for reactions that introduce bromine atoms
            if not is_bromination:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has more bromine atoms than reactants
                reactants_br_count = sum(reactant.count("Br") for reactant in reactants)
                product_br_count = product.count("Br")

                if product_br_count > reactants_br_count:
                    is_bromination = True
                    print(f"Bromination detected by atom count: {rsmi}")

                # Also check if a bromine atom is introduced at a specific position
                # This is important for strategic bromination
                if not is_bromination and "Br" in product:
                    # If the product has a bromine and it's part of a synthetic strategy
                    # where the bromine is retained through multiple steps
                    if aromatic_bromide_in_product(route):
                        is_bromination = True
                        print(f"Strategic bromination detected by product analysis: {rsmi}")

            if is_bromination:
                print(f"Strategic bromination step detected: {rsmi}")
                bromination_found = True

        # Recursively check children
        for child in node.get("children", []):
            dfs(child)

    # Start DFS from the root
    dfs(route)

    return bromination_found
