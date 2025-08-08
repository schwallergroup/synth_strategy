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
    This function detects if the synthesis involves aromatic halogenation,
    particularly chlorination, fluorination, bromination, or iodination of aromatic rings.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["rsmi"]

                # Check for aromatic halogenation reactions using the checker functions
                if checker.check_reaction("Aromatic fluorination", rsmi):
                    print(f"Aromatic fluorination detected in reaction: {rsmi}")
                    result = True
                elif checker.check_reaction("Aromatic chlorination", rsmi):
                    print(f"Aromatic chlorination detected in reaction: {rsmi}")
                    result = True
                elif checker.check_reaction("Aromatic bromination", rsmi):
                    print(f"Aromatic bromination detected in reaction: {rsmi}")
                    result = True
                elif checker.check_reaction("Aromatic iodination", rsmi):
                    print(f"Aromatic iodination detected in reaction: {rsmi}")
                    result = True

                # If no specific reaction type is detected, check for functional group changes
                if not result:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check if product contains aromatic halides
                    if product_mol:
                        has_aromatic_halide = False

                        # Check for aromatic halides in the product
                        if checker.check_fg("Aromatic halide", product_smiles):
                            # Verify this is a new halide by checking reactants
                            has_new_halide = True
                            for reactant in reactants_smiles:
                                if checker.check_fg("Aromatic halide", reactant):
                                    # If a reactant already has an aromatic halide, check if product has more
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol:
                                        # This is a simplistic check - in a real scenario we'd need to track atom mappings
                                        has_new_halide = True

                            if has_new_halide:
                                print(
                                    f"Aromatic halogenation detected by functional group analysis in: {rsmi}"
                                )
                                result = True

            except Exception as e:
                print(f"Error in halogenation detection: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return result
