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
    This function detects synthesis routes involving trifluoromethyl-containing compounds.
    It checks both molecule nodes and reaction nodes for the presence of trifluoromethyl groups.
    """
    has_trifluoromethyl = False

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoromethyl

        # Process molecule nodes
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]
                if checker.check_fg("Trifluoro group", mol_smiles):
                    print(f"Detected trifluoromethyl group in molecule: {mol_smiles}")
                    has_trifluoromethyl = True
            except Exception as e:
                print(f"Error processing molecule node: {e}")

        # Process reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant contains trifluoromethyl group
                for reactant in reactants:
                    if checker.check_fg("Trifluoro group", reactant):
                        print(f"Detected trifluoromethyl group in reactant: {reactant}")
                        has_trifluoromethyl = True

                # Check if product contains trifluoromethyl group
                if checker.check_fg("Trifluoro group", product):
                    print(f"Detected trifluoromethyl group in product: {product}")
                    has_trifluoromethyl = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_trifluoromethyl
