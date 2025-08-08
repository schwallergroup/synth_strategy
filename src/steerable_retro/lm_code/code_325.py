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
    Detects if the synthetic route involves a nitro to amine conversion.
    """

    def dfs_traverse(node, depth=0):
        # Check reaction nodes
        if node["type"] == "reaction":
            try:
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    # Check if this is a nitro reduction reaction
                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        print(f"Found nitro reduction reaction at depth {depth}: {rsmi}")
                        return True

                    # Alternative check: look at reactants and products
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if reactant has nitro group and product has amine group
                    reactant_has_nitro = any(checker.check_fg("Nitro group", r) for r in reactants)
                    product_has_amine = checker.check_fg("Primary amine", product)

                    if reactant_has_nitro and product_has_amine:
                        print(f"Found nitro to amine conversion at depth {depth}: {rsmi}")
                        return True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Recursively check children
        for child in node.get("children", []):
            if dfs_traverse(child, depth + 1):
                return True

        return False

    # Start traversal
    return dfs_traverse(route)
