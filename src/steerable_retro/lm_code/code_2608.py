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
    This function detects a synthetic strategy involving thiazole heterocycle formation
    from fragment coupling.
    """
    found_thiazole_formation = False

    def dfs_traverse(node):
        nonlocal found_thiazole_formation

        if node["type"] == "reaction":
            metadata = node.get("metadata", {})
            rsmi = metadata.get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants_smiles = parts[0].split(".")
            product_smiles = parts[2]

            # Check if this is a thiazole formation reaction
            if checker.check_reaction("thiazole", rsmi):
                print(f"Found thiazole formation reaction: {rsmi}")
                found_thiazole_formation = True
                return

            # Alternative check: look for thiazole ring formation
            product_has_thiazole = checker.check_ring("thiazole", product_smiles)

            if product_has_thiazole:
                # Check if any reactant already has a thiazole ring
                reactants_have_thiazole = any(
                    checker.check_ring("thiazole", r) for r in reactants_smiles if r
                )

                if not reactants_have_thiazole:
                    print(f"Found thiazole formation from fragment coupling")
                    print(f"Reactants: {reactants_smiles}")
                    print(f"Product: {product_smiles}")
                    found_thiazole_formation = True

                    # Additional check for specific reaction patterns
                    # Check for common thiazole formation reactions if not already detected
                    if (
                        checker.check_reaction("benzothiazole", rsmi)
                        or checker.check_reaction(
                            "benzothiazole_derivatives_carboxylic-acid/ester", rsmi
                        )
                        or checker.check_reaction("benzothiazole_derivatives_aldehyde", rsmi)
                    ):
                        print(f"Confirmed specific thiazole formation reaction type")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_thiazole_formation
