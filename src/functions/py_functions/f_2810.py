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


def main(route):
    """
    Detects if the synthetic route follows a linear synthesis strategy
    (each step has only one non-commercial reactant)

    A linear synthesis strategy means that at each reaction step, only one
    reactant is a complex intermediate (non-commercial), while other reactants
    are commercial reagents. This creates a linear sequence of reactions without
    convergent steps.

    Args:
        route: A dictionary representing the synthesis route

    Returns:
        bool: True if the route follows a linear synthesis strategy, False otherwise
    """
    is_linear = True
    reaction_count = 0

    # Dictionary to store commercial status of molecules by SMILES
    commercial_status = {}

    # First pass: identify all commercial molecules
    def identify_commercial(node):
        if node["type"] == "mol":
            smiles = node["smiles"]
            if "in_stock" in node and node["in_stock"]:
                commercial_status[smiles] = True
            else:
                commercial_status[smiles] = False

        for child in node.get("children", []):
            identify_commercial(child)

    # Check if a reactant is commercial based on its SMILES
    def is_commercial_reactant(reactant_smiles, reaction_node):
        # Small molecules are typically reagents
        mol = Chem.MolFromSmiles(reactant_smiles)
        if mol is None:
            print(f"Warning: Could not parse SMILES: {reactant_smiles}")
            return True  # Assume commercial if can't parse

        if mol.GetNumAtoms() <= 3:  # Very small molecules are likely commercial
            return True

        # Check if we've already determined commercial status
        if reactant_smiles in commercial_status:
            return commercial_status[reactant_smiles]

        # Check if this reactant appears as a child node of the reaction
        for child in reaction_node.get("children", []):
            if child["type"] == "mol" and child["smiles"] == reactant_smiles:
                if "in_stock" in child and child["in_stock"]:
                    commercial_status[reactant_smiles] = True
                    return True
                else:
                    commercial_status[reactant_smiles] = False
                    return False

        # If we can't find it in children, check for atom-mapped version
        try:
            # Remove atom mapping for comparison
            unmapped_smiles = re.sub(r":[0-9]+", "", reactant_smiles)
            for child in reaction_node.get("children", []):
                if child["type"] == "mol":
                    child_unmapped = re.sub(r":[0-9]+", "", child["smiles"])
                    if child_unmapped == unmapped_smiles:
                        if "in_stock" in child and child["in_stock"]:
                            commercial_status[reactant_smiles] = True
                            return True
                        else:
                            commercial_status[reactant_smiles] = False
                            return False
        except Exception as e:
            print(f"Error comparing atom-mapped SMILES: {e}")

        # Default to assuming it's a commercial reagent if we can't determine
        commercial_status[reactant_smiles] = True
        return True

    # Second pass: check if synthesis is linear
    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            reaction_count += 1
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]

            # Count non-commercial reactants
            reactants = reactants_part.split(".")
            non_commercial_count = 0

            for reactant in reactants:
                if not is_commercial_reactant(reactant, node):
                    non_commercial_count += 1

            # If more than one non-commercial reactant, it's not a linear synthesis
            if non_commercial_count > 1:
                print(
                    f"Found convergent step with {non_commercial_count} non-commercial reactants"
                )
                is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child)

    # First identify all commercial molecules
    identify_commercial(route)

    # Then check if synthesis is linear
    dfs_traverse(route)

    # Need at least 2 reactions to determine if it's linear
    if reaction_count < 2:
        print("Route has fewer than 2 reactions, not considering it linear")
        return False

    if is_linear:
        print(
            f"Route follows a linear synthesis strategy with {reaction_count} reactions"
        )

    return is_linear
