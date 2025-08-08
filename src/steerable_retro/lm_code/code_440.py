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
    This function detects if the synthesis involves borylation of an aryl halide
    as preparation for a coupling reaction.
    """
    has_borylation = False

    def dfs_traverse(node):
        nonlocal has_borylation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a borylation reaction using the reaction checker
            is_borylation_reaction = (
                checker.check_reaction("Preparation of boronic acids", rsmi)
                or checker.check_reaction(
                    "Preparation of boronic acids without boronic ether", rsmi
                )
                or checker.check_reaction("Preparation of boronic ethers", rsmi)
                or checker.check_reaction(
                    "Preparation of boronic acids from trifluoroborates", rsmi
                )
            )

            if is_borylation_reaction:
                print(f"Found borylation reaction by reaction type: {rsmi}")
                has_borylation = True
            else:
                # Check for aryl halide in reactants
                has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

                # Check for boronic acid/ester in product
                has_boronic_product = checker.check_fg("Boronic acid", product) or checker.check_fg(
                    "Boronic ester", product
                )

                if has_aryl_halide and has_boronic_product:
                    print(f"Found borylation reaction by functional group analysis: {rsmi}")
                    has_borylation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_borylation
