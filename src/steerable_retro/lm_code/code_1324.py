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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    Detects if the final step (depth 0) in the synthesis is an ester hydrolysis.
    """
    final_step_is_hydrolysis = False

    def is_final_step(node, route):
        """Helper function to determine if this is the final step when depth isn't available"""
        try:
            if node["type"] == "reaction":
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check if this product matches the target molecule
                from rdkit import Chem

                product_mol = Chem.MolFromSmiles(product)
                route_mol = Chem.MolFromSmiles(route["smiles"])
                if product_mol and route_mol:
                    return Chem.MolToSmiles(product_mol) == Chem.MolToSmiles(route_mol)
        except Exception as e:
            print(f"Error in is_final_step: {e}")
        return False

    def dfs_traverse(node, depth=None, is_root_child=False):
        nonlocal final_step_is_hydrolysis

        if node["type"] == "reaction":
            # Get depth from metadata if available
            if depth is None:
                depth = node.get("metadata", {}).get("depth")

            # Check if it's the final step (depth 0 or direct child of root or produces final product)
            is_final = (depth == 0) or is_root_child or is_final_step(node, route)

            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}, is_final={is_final}, rsmi: {rsmi}")

                if is_final:
                    # Check if any reactant contains an ester
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                    # Check if product contains carboxylic acid
                    has_acid = checker.check_fg("Carboxylic acid", product)

                    # Verify this is a hydrolysis reaction
                    is_hydrolysis = (
                        checker.check_reaction(
                            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                        )
                        or checker.check_reaction(
                            "Ester saponification (methyl deprotection)", rsmi
                        )
                        or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    )

                    print(
                        f"Final step check: has_ester={has_ester}, has_acid={has_acid}, is_hydrolysis={is_hydrolysis}"
                    )

                    if has_ester and has_acid and is_hydrolysis:
                        print("Found ester hydrolysis as final step")
                        final_step_is_hydrolysis = True
                    elif has_ester and has_acid:
                        # Additional check for ester to acid conversion even if reaction type check failed
                        print("Found ester to acid conversion, checking manually")
                        # Look for ester in reactants and corresponding acid in product using atom mapping
                        from rdkit import Chem

                        for reactant in reactants:
                            if checker.check_fg("Ester", reactant):
                                r_mol = Chem.MolFromSmiles(reactant)
                                p_mol = Chem.MolFromSmiles(product)
                                if r_mol and p_mol:
                                    # If we have an ester to acid conversion, it's likely a hydrolysis
                                    final_step_is_hydrolysis = True
                                    print("Confirmed ester hydrolysis through manual check")
                                    break
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1 if depth is not None else None)

    # First check direct children of the root (final steps)
    for child in route.get("children", []):
        dfs_traverse(child, depth=0, is_root_child=True)

    # If no hydrolysis found in direct children, do a full traversal
    if not final_step_is_hydrolysis:
        dfs_traverse(route)

    return final_step_is_hydrolysis
