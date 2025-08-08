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
    This function detects if the synthetic route involves a sequence of
    nitro group introduction followed by reduction to amine.
    """
    nitration_steps = []
    reduction_steps = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitration reactions
                if (
                    checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                    or checker.check_reaction("Non-aromatic nitration with HNO3", rsmi)
                ):
                    nitration_steps.append(depth)
                    print(f"Detected nitration at depth {depth}")

                # If not a specific nitration reaction, check for nitro group appearance
                elif checker.check_fg("Nitro group", product_smiles):
                    # Check if nitro group wasn't in reactants
                    nitro_in_reactants = False
                    for reactant_smiles in reactants_smiles:
                        if checker.check_fg("Nitro group", reactant_smiles):
                            nitro_in_reactants = True
                            break

                    if not nitro_in_reactants:
                        nitration_steps.append(depth)
                        print(f"Detected nitro group introduction at depth {depth}")

                # Check for nitro reduction to amine
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    reduction_steps.append(depth)
                    print(f"Detected nitro reduction at depth {depth}")
                else:
                    # Alternative check: amine appears in product while nitro was in reactants
                    if checker.check_fg("Primary amine", product_smiles):
                        # Check if amine wasn't in reactants but nitro was
                        amine_in_reactants = False
                        nitro_in_reactants = False

                        for reactant_smiles in reactants_smiles:
                            if checker.check_fg("Primary amine", reactant_smiles):
                                amine_in_reactants = True
                            if checker.check_fg("Nitro group", reactant_smiles):
                                nitro_in_reactants = True

                        if not amine_in_reactants and nitro_in_reactants:
                            reduction_steps.append(depth)
                            print(f"Detected nitro reduction at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if there's a nitration followed by reduction
    for nitration_depth in nitration_steps:
        for reduction_depth in reduction_steps:
            # In retrosynthetic traversal, lower depth = later in synthesis (earlier in retrosynthesis)
            if (
                reduction_depth < nitration_depth
            ):  # Check if reduction happens after nitration in forward synthesis
                print("Detected nitro-reduction-amination sequence")
                return True

    return False
