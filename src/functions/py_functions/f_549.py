#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis includes anhydride chemistry.
    """
    has_anhydride = False

    def dfs_traverse(node):
        nonlocal has_anhydride

        # Check molecule nodes for anhydride functional group
        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]
            if checker.check_fg("Anhydride", mol_smiles):
                print(f"Anhydride functional group detected in molecule: {mol_smiles}")
                has_anhydride = True

        # Check reaction nodes for anhydride chemistry
        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for anhydride pattern in reactants
            for reactant in reactants:
                if reactant and checker.check_fg("Anhydride", reactant):
                    print(
                        f"Anhydride functional group detected in reactant: {reactant}"
                    )
                    has_anhydride = True

            # Check for anhydride pattern in product
            if product and checker.check_fg("Anhydride", product):
                print(f"Anhydride functional group detected in product: {product}")
                has_anhydride = True

            # Check for specific anhydride-related reactions
            if any(
                checker.check_reaction(rxn_name, rsmi)
                for rxn_name in [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acylation of secondary amines with anhydrides",
                ]
            ):
                print(f"Anhydride-related reaction detected")
                has_anhydride = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_anhydride
