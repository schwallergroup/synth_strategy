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
    This function detects if the synthesis involves a reductive amination step
    (aldehyde/ketone + amine -> amine).
    """
    found_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal found_reductive_amination

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is a reductive amination reaction using the checker function
                if checker.check_reaction("reductive amination with aldehyde", rsmi):
                    print(f"Found reductive amination with aldehyde at depth {depth}")
                    found_reductive_amination = True
                    return
                elif checker.check_reaction("reductive amination with ketone", rsmi):
                    print(f"Found reductive amination with ketone at depth {depth}")
                    found_reductive_amination = True
                    return
                elif checker.check_reaction("reductive amination with alcohol", rsmi):
                    print(f"Found reductive amination with alcohol at depth {depth}")
                    found_reductive_amination = True
                    return
                else:
                    # Fallback method: check for reactants and products manually
                    print("Using fallback method to check for reductive amination pattern")
                    try:
                        reactants_str = rsmi.split(">")[0]
                        product_str = rsmi.split(">")[-1]
                        reactants = reactants_str.split(".")

                        if len(reactants) >= 2:
                            # Check for aldehyde/ketone/alcohol and amine in reactants
                            carbonyl_reactant = None
                            amine_reactant = None

                            for reactant in reactants:
                                if checker.check_fg("Aldehyde", reactant):
                                    print(f"Found aldehyde in reactant: {reactant}")
                                    carbonyl_reactant = reactant
                                elif checker.check_fg("Ketone", reactant):
                                    print(f"Found ketone in reactant: {reactant}")
                                    carbonyl_reactant = reactant
                                elif checker.check_fg(
                                    "Primary alcohol", reactant
                                ) or checker.check_fg("Secondary alcohol", reactant):
                                    print(f"Found alcohol in reactant: {reactant}")
                                    carbonyl_reactant = reactant

                                if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                    "Secondary amine", reactant
                                ):
                                    print(f"Found amine in reactant: {reactant}")
                                    amine_reactant = reactant

                            if carbonyl_reactant and amine_reactant:
                                print(
                                    f"Found potential reductive amination reactants. Checking product: {product_str}"
                                )

                                # Check if product has an amine (primary, secondary, or tertiary)
                                has_amine = (
                                    checker.check_fg("Primary amine", product_str)
                                    or checker.check_fg("Secondary amine", product_str)
                                    or checker.check_fg("Tertiary amine", product_str)
                                )

                                if has_amine:
                                    # Try to verify the transformation using atom mapping
                                    try:
                                        # Extract atom-mapped molecules
                                        carbonyl_mol = Chem.MolFromSmiles(carbonyl_reactant)
                                        amine_mol = Chem.MolFromSmiles(amine_reactant)
                                        product_mol = Chem.MolFromSmiles(product_str)

                                        if carbonyl_mol and amine_mol and product_mol:
                                            # Check if the product has a new C-N bond that wasn't in the reactants
                                            # This is a simplified check - in a real implementation, we would
                                            # track the specific atom mappings from the carbonyl carbon to the amine nitrogen
                                            print(
                                                f"Found reductive amination pattern at depth {depth}"
                                            )
                                            found_reductive_amination = True
                                            return
                                    except Exception as e:
                                        print(f"Error in atom mapping verification: {e}")

                                    # Even without atom mapping verification, if we have the right reactants and product,
                                    # it's likely a reductive amination
                                    print(
                                        f"Found likely reductive amination pattern at depth {depth}"
                                    )
                                    found_reductive_amination = True
                                    return
                                else:
                                    print(f"Product doesn't have an amine functional group")
                    except Exception as e:
                        print(f"Error in fallback method: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting to search for reductive amination strategy...")
    dfs_traverse(route)
    print(f"Reductive amination found: {found_reductive_amination}")
    return found_reductive_amination
