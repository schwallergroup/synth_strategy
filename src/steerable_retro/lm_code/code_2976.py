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
    Detects if cyano groups are retained throughout the synthesis.

    This function checks if cyano groups (nitrile groups) are preserved in the
    synthetic route. It verifies that once a cyano group appears in a molecule,
    it is retained in all subsequent steps of the synthesis.
    """
    # Track molecules that should have cyano groups
    molecules_with_cyano = set()
    molecules_without_cyano = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            smiles = node["smiles"]

            # Check if this molecule has a cyano group
            has_cyano = checker.check_fg("Nitrile", smiles)

            # If it's not a starting material, check if it should have a cyano group
            if not node.get("in_stock", False):
                if smiles in molecules_with_cyano and not has_cyano:
                    print(f"Molecule lost cyano group: {smiles}")
                    molecules_without_cyano.add(smiles)
                elif has_cyano:
                    molecules_with_cyano.add(smiles)
            # If it's a starting material with a cyano group, track it
            elif has_cyano:
                molecules_with_cyano.add(smiles)

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if any reactant has a cyano group
            reactants = reactants_part.split(".")
            for reactant in reactants:
                if checker.check_fg("Nitrile", reactant):
                    # If a reactant has cyano, the product should also have it
                    if not checker.check_fg("Nitrile", product_part):
                        print(f"Cyano group lost in reaction: {rsmi}")
                        molecules_without_cyano.add(product_part)
                    else:
                        molecules_with_cyano.add(product_part)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # If any molecule that should have a cyano group doesn't have one, return False
    return len(molecules_without_cyano) == 0
