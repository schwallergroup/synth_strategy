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
    Detects if an indazole core is maintained throughout the synthesis.

    This function checks if once an indazole core appears in the synthesis route,
    it is maintained in all subsequent steps toward the final product.
    """
    # Check if the target molecule contains indazole
    if route["type"] == "mol" and "smiles" in route:
        target_has_indazole = checker.check_ring("indazole", route["smiles"])
        if not target_has_indazole:
            print(f"Target molecule does not contain indazole: {route['smiles']}")
            return False
        print(f"Target molecule contains indazole: {route['smiles']}")
    else:
        print("Route does not start with a molecule node")
        return False

    # Track the synthetic pathway
    def check_pathway(node, indazole_required=True):
        # For molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_indazole = checker.check_ring("indazole", mol_smiles)

            # If this is a starting material
            if node.get("in_stock", False):
                # Starting materials don't need to have indazole
                print(f"Starting material: {mol_smiles}, has indazole: {has_indazole}")
                return True

            # If indazole is required but not present
            if indazole_required and not has_indazole:
                print(f"Indazole core not maintained in: {mol_smiles}")
                return False

            # If no children, this is a leaf node (starting material)
            if not node.get("children", []):
                return True

            # Check all children (reactions)
            for child in node.get("children", []):
                if not check_pathway(child, has_indazole):
                    return False
            return True

        # For reaction nodes
        elif node["type"] == "reaction":
            # Get reaction SMILES if available
            rxn_smiles = node.get("metadata", {}).get("rsmi", "")
            print(f"Checking reaction: {rxn_smiles}")

            # For each reactant (child molecule), at least one should maintain indazole if required
            reactants_with_indazole = 0
            total_reactants = 0

            for child in node.get("children", []):
                if child["type"] == "mol":
                    total_reactants += 1
                    if checker.check_ring("indazole", child["smiles"]):
                        reactants_with_indazole += 1
                        # Continue checking this branch
                        if not check_pathway(child, True):
                            return False
                    else:
                        # This reactant doesn't have indazole, but that's okay if it's not the main one
                        if not check_pathway(child, False):
                            return False

            # If indazole is required in this reaction, at least one reactant should have it
            if indazole_required and reactants_with_indazole == 0 and total_reactants > 0:
                print(f"No reactants contain indazole in reaction: {rxn_smiles}")
                return False

            return True

        return True

    # Start checking from the root
    return check_pathway(route, True)
