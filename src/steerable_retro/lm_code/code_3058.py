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
    """Check if the route involves formation of a tertiary alcohol."""
    tertiary_alcohol_found = False

    def dfs(node, depth=0):
        nonlocal tertiary_alcohol_found

        if tertiary_alcohol_found:
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check if the product contains a tertiary alcohol
            if checker.check_fg("Tertiary alcohol", product):
                # Check reactants to confirm this is a formation reaction
                reactants = rsmi.split(">")[0].split(".")
                reactant_has_tertiary_alcohol = any(
                    checker.check_fg("Tertiary alcohol", r) for r in reactants
                )

                if not reactant_has_tertiary_alcohol:
                    tertiary_alcohol_found = True
                    print(f"Found tertiary alcohol formation reaction: {rsmi}")
                    return

                # Also check for specific reactions that form tertiary alcohols
                if checker.check_reaction("Grignard from ketone to alcohol", rsmi):
                    tertiary_alcohol_found = True
                    print(f"Found tertiary alcohol formation via Grignard: {rsmi}")
                    return

        # Also check the target molecule for tertiary alcohols
        if node["type"] == "mol" and node.get("children", []) and not node.get("in_stock", False):
            mol_smiles = node["smiles"]
            if checker.check_fg("Tertiary alcohol", mol_smiles):
                tertiary_alcohol_found = True
                print(f"Found tertiary alcohol in non-starting material: {mol_smiles}")
                return

        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return tertiary_alcohol_found
