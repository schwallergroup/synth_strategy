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


def main(route):
    """
    Detects if the synthesis follows a linear strategy (as opposed to convergent).

    A linear synthesis strategy is one where each reaction step has at most one
    non-commercial (not in_stock) reactant. If any reaction has more than one
    non-commercial reactant, the synthesis is considered convergent.

    Args:
        route: A dictionary representing the synthesis route following the SynthesisRoute schema

    Returns:
        bool: True if the synthesis follows a linear strategy, False otherwise
    """
    # In a linear synthesis, each reaction typically has only one non-commercial reactant
    linear_strategy = True

    def dfs_traverse(node, depth=0):
        nonlocal linear_strategy

        # If we've already determined it's not linear, no need to continue
        if not linear_strategy:
            return

        if not isinstance(node, dict) or "type" not in node:
            print(f"Invalid node structure at depth {depth}")
            return

        if node["type"] == "reaction":
            try:
                # Get the reaction SMILES from metadata
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    print(f"Warning: Missing reaction SMILES at depth {depth}")
                    return

                # Extract reactants from the reaction SMILES
                reactants_part = rsmi.split(">")[0]
                reactants_smiles = reactants_part.split(".")

                # Count non-commercial reactants
                non_commercial_reactants = 0
                non_commercial_smiles = []

                # Process each child node (reactants)
                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        # This is a non-commercial reactant
                        child_smiles = child.get("smiles", "")
                        child_mol = Chem.MolFromSmiles(child_smiles)

                        if child_mol:
                            # Check if this child corresponds to any reactant in the reaction
                            is_reactant = False
                            for reactant_smiles in reactants_smiles:
                                # Remove atom mapping for comparison
                                clean_reactant_smiles = reactant_smiles
                                for atom in clean_reactant_smiles:
                                    if ":" in clean_reactant_smiles:
                                        clean_reactant_smiles = clean_reactant_smiles.replace(
                                            atom, atom.split(":")[0]
                                        )

                                reactant_mol = Chem.MolFromSmiles(clean_reactant_smiles)
                                if reactant_mol:
                                    # Use isomorphism to check if molecules are the same
                                    if Chem.MolToSmiles(child_mol) == Chem.MolToSmiles(
                                        reactant_mol
                                    ):
                                        is_reactant = True
                                        break

                            if is_reactant:
                                non_commercial_reactants += 1
                                non_commercial_smiles.append(child_smiles)

                # If more than one non-commercial reactant, it's convergent
                if non_commercial_reactants > 1:
                    linear_strategy = False
                    print(
                        f"Found convergent step at depth {depth} with {non_commercial_reactants} non-commercial reactants"
                    )
                    print(f"Reaction: {rsmi}")
                    print(f"Non-commercial reactants: {', '.join(non_commercial_smiles)}")
                    return  # Early return once we find a convergent step

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {str(e)}")
                # Continue processing other nodes even if this one fails

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return linear_strategy
