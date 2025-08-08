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
    Detects if the hydroxymethyl group (CH2OH) is preserved throughout the synthesis.
    A preservation strategy means the hydroxymethyl group is present in the final product
    and is maintained intact through the synthesis route.
    """
    # Check if the final product (root node) has a hydroxymethyl group
    if route["type"] != "mol" or "smiles" not in route:
        print("Root node is not a molecule or missing SMILES")
        return False

    target_mol_smiles = route["smiles"]
    target_mol = Chem.MolFromSmiles(target_mol_smiles)

    # Check if hydroxymethyl exists in the final product
    has_hydroxymethyl = checker.check_fg("Primary alcohol", target_mol_smiles)
    if not has_hydroxymethyl:
        print("No hydroxymethyl group in the final product")
        return False

    # Get the atom indices of the hydroxymethyl group in the final product
    hydroxymethyl_indices = checker.get_fg_atom_indices("Primary alcohol", target_mol_smiles)
    if not hydroxymethyl_indices:
        print("Could not identify hydroxymethyl atom indices")
        return False

    # Track if the hydroxymethyl group is preserved throughout the synthesis
    preserved = [True]  # Using list to allow modification in nested function

    def dfs_traverse(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this reaction preserves the hydroxymethyl group
            # by ensuring it's not directly involved in the reaction
            if not checker.check_fg("Primary alcohol", rsmi.split(">")[-1]):
                print(f"Hydroxymethyl not preserved in reaction: {rsmi}")
                preserved[0] = False
                return

            # Continue checking children (reactants)
            for child in node.get("children", []):
                dfs_traverse(child)

        elif node["type"] == "mol" and node != route:  # Skip the target molecule (already checked)
            if "smiles" in node and not node.get("in_stock", False):
                mol_smiles = node["smiles"]
                # Check if intermediate also has the hydroxymethyl group
                if not checker.check_fg("Primary alcohol", mol_smiles):
                    print(f"Hydroxymethyl not found in intermediate: {mol_smiles}")
                    preserved[0] = False

    # Start traversal from the root
    for child in route.get("children", []):
        dfs_traverse(child)

    if preserved[0]:
        print("Hydroxymethyl preserved throughout synthesis")
        return True
    return False
