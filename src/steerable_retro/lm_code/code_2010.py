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
    This function detects a strategy where a tosylate is used as a leaving group
    for displacement by an azide nucleophile.
    """
    # Initialize flag
    has_tosylate_azide_displacement = False

    def dfs_traverse(node, depth=0):
        nonlocal has_tosylate_azide_displacement

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check if this is a nucleophilic substitution reaction that could involve azide formation
                    if (
                        checker.check_reaction("Alcohol to azide", rsmi)
                        or checker.check_reaction("Formation of Azides from halogens", rsmi)
                        or checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                    ):

                        # Check if one reactant has tosylate
                        has_tosylate = False
                        for reactant in reactants:
                            if reactant.strip():
                                if checker.check_fg("Tosylate", reactant):
                                    has_tosylate = True
                                    print(f"Found tosylate in reactant: {reactant}")
                                    break

                        # Check if product has azide
                        has_azide = checker.check_fg("Azide", product)
                        if has_azide:
                            print(f"Found azide in product: {product}")

                        # If both conditions are met, we've found the strategy
                        if has_tosylate and has_azide:
                            has_tosylate_azide_displacement = True
                            print(f"Detected tosylate displacement by azide at depth {depth}")

                    # Also check for any reaction where tosylate is converted to azide
                    # This is a more general check that doesn't rely on specific reaction types
                    else:
                        has_tosylate_in_reactants = False
                        has_azide_in_product = checker.check_fg("Azide", product)

                        for reactant in reactants:
                            if reactant.strip() and checker.check_fg("Tosylate", reactant):
                                has_tosylate_in_reactants = True
                                break

                        if has_tosylate_in_reactants and has_azide_in_product:
                            # Check if the product doesn't have tosylate, indicating it was displaced
                            if not checker.check_fg("Tosylate", product):
                                has_tosylate_azide_displacement = True
                                print(
                                    f"Detected general tosylate to azide transformation at depth {depth}"
                                )
                except Exception as e:
                    print(f"Error in reaction analysis: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    if has_tosylate_azide_displacement:
        print("Detected tosylate-azide displacement strategy")

    return has_tosylate_azide_displacement
