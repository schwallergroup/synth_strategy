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
    Detects the use of Fischer indole synthesis to construct the indole core.
    """
    fischer_indole_found = False

    def dfs_traverse(node):
        nonlocal fischer_indole_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # First check if this is a Fischer indole reaction
                if checker.check_reaction("Fischer indole", rsmi):
                    print(f"Detected Fischer indole reaction directly: {rsmi}")
                    fischer_indole_found = True
                else:
                    # Check for required components of Fischer indole synthesis
                    reactant_has_hydrazine = False
                    reactant_has_carbonyl = False

                    for reactant in reactants:
                        # Check for hydrazine/hydrazone functional group
                        if checker.check_fg("Hydrazine", reactant) or checker.check_fg(
                            "Hydrazone", reactant
                        ):
                            print(f"Found hydrazine/hydrazone reactant: {reactant}")
                            reactant_has_hydrazine = True

                        # Check for ketone or aldehyde (carbonyl compounds)
                        if checker.check_fg("Ketone", reactant) or checker.check_fg(
                            "Aldehyde", reactant
                        ):
                            print(f"Found carbonyl reactant: {reactant}")
                            reactant_has_carbonyl = True

                    # Check for indole ring in product
                    product_has_indole = checker.check_ring("indole", product)
                    if product_has_indole:
                        print(f"Found indole in product: {product}")

                    # If we have all components of Fischer indole synthesis
                    if reactant_has_hydrazine and reactant_has_carbonyl and product_has_indole:
                        print(f"Detected Fischer indole synthesis pattern in reaction: {rsmi}")
                        fischer_indole_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return fischer_indole_found
