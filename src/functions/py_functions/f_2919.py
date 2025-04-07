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
    Detects if the synthetic route follows a linear synthesis strategy (no convergent steps).

    A linear synthesis strategy has only one non-commercial (synthesized) reactant per step.
    If any reaction step has multiple synthesized reactants, it's considered convergent.
    """
    is_linear = True

    # First, collect all molecule nodes in the route for quick lookup
    molecule_nodes = {}

    def collect_molecules(node):
        if node["type"] == "mol" and "smiles" in node:
            # Canonicalize SMILES for consistent comparison
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    canonical_smiles = Chem.MolToSmiles(mol)
                    molecule_nodes[canonical_smiles] = node
            except Exception as e:
                print(f"Error processing molecule SMILES: {e}")

        for child in node.get("children", []):
            collect_molecules(child)

    collect_molecules(route)

    def dfs_traverse(node):
        nonlocal is_linear

        if not is_linear:  # Early return if we already found a convergent step
            return

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Count synthesized reactants (those that are in our synthesis tree and not in_stock)
            synthesized_reactants = 0

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        canonical_smiles = Chem.MolToSmiles(mol)

                        # Check if this reactant is in our synthesis tree
                        if canonical_smiles in molecule_nodes:
                            reactant_node = molecule_nodes[canonical_smiles]

                            # If it's not marked as in_stock, it's a synthesized reactant
                            if not reactant_node.get("in_stock", False):
                                synthesized_reactants += 1
                                print(f"Found synthesized reactant: {canonical_smiles}")
                except Exception as e:
                    print(f"Error processing reactant SMILES: {e}")

            # If more than one synthesized reactant, it's a convergent step
            if synthesized_reactants > 1:
                print(
                    f"Convergent step detected with {synthesized_reactants} synthesized reactants"
                )
                is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
