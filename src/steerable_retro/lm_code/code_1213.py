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
    Detects if a bromine substituent is maintained throughout the entire synthesis.

    This function checks if bromine atoms are present in the main synthetic pathway,
    excluding starting materials that don't contribute bromine to the synthesis.
    """
    # Track if we've found at least one molecule with bromine
    found_bromine = False
    # Track if bromine is maintained throughout synthesis
    maintains_bromine = True

    def dfs_traverse(node, depth=0):
        nonlocal found_bromine, maintains_bromine

        # Check molecule nodes
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            has_bromine = (
                checker.check_fg("Aromatic halide", mol_smiles)
                or checker.check_fg("Primary halide", mol_smiles)
                or checker.check_fg("Secondary halide", mol_smiles)
                or checker.check_fg("Tertiary halide", mol_smiles)
                or checker.check_fg("Alkenyl halide", mol_smiles)
            )

            # If this is the target molecule (depth 0), it must have bromine
            if depth == 0:
                if not has_bromine:
                    print(f"Target molecule doesn't have bromine: {mol_smiles}")
                    maintains_bromine = False
                else:
                    found_bromine = True

            # For intermediate molecules (not starting materials), check bromine retention
            elif not node.get("in_stock", False) and found_bromine:
                if not has_bromine:
                    print(f"Intermediate molecule lost bromine: {mol_smiles}")
                    maintains_bromine = False

        # Check reaction nodes to ensure bromine transfers from reactants to products
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has bromine
            reactant_has_bromine = any(
                checker.check_fg("Aromatic halide", r)
                or checker.check_fg("Primary halide", r)
                or checker.check_fg("Secondary halide", r)
                or checker.check_fg("Tertiary halide", r)
                or checker.check_fg("Alkenyl halide", r)
                for r in reactants
            )

            # Check if product has bromine
            product_has_bromine = (
                checker.check_fg("Aromatic halide", product)
                or checker.check_fg("Primary halide", product)
                or checker.check_fg("Secondary halide", product)
                or checker.check_fg("Tertiary halide", product)
                or checker.check_fg("Alkenyl halide", product)
            )

            # If reactants have bromine but product doesn't, bromine was lost
            if reactant_has_bromine and not product_has_bromine:
                print(f"Bromine lost in reaction: {rsmi}")
                maintains_bromine = False

            # If reactants have bromine, mark that we've found bromine in the synthesis
            if reactant_has_bromine:
                found_bromine = True

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # We need to have found at least one bromine and maintained it throughout
    return found_bromine and maintains_bromine
