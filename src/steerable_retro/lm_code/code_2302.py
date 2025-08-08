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
    This function detects a synthetic strategy involving multiple Boc deprotection steps.
    """
    boc_deprotection_count = 0

    def dfs_traverse(node):
        nonlocal boc_deprotection_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a Boc deprotection reaction using reaction type
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    boc_deprotection_count += 1
                    print(f"Found Boc deprotection step: {rsmi}")
                # If reaction type check fails, try checking for Boc group disappearance
                else:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Count actual Boc groups in reactants and product
                    boc_count_in_reactants = 0
                    for r in reactants:
                        try:
                            # Get all instances of Boc groups in each reactant
                            boc_indices = checker.get_fg_atom_indices("Boc", r)
                            if boc_indices:
                                boc_count_in_reactants += len(boc_indices)
                        except:
                            # Handle potential errors in functional group detection
                            pass

                    boc_count_in_product = 0
                    try:
                        # Get all instances of Boc groups in the product
                        boc_indices = checker.get_fg_atom_indices("Boc", product)
                        if boc_indices:
                            boc_count_in_product = len(boc_indices)
                    except:
                        # Handle potential errors in functional group detection
                        pass

                    # Calculate how many Boc groups were removed
                    boc_groups_removed = boc_count_in_reactants - boc_count_in_product

                    if boc_groups_removed > 0:
                        # Add the number of Boc groups removed
                        boc_deprotection_count += boc_groups_removed
                        print(f"Found Boc deprotection step (detected by FG change): {rsmi}")
                        print(f"Number of Boc groups removed: {boc_groups_removed}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total Boc deprotection steps: {boc_deprotection_count}")
    # Strategy requires at least 2 Boc deprotection steps
    return boc_deprotection_count >= 2
