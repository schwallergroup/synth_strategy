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
    Detects if the synthesis uses a convergent approach by checking if
    the final step combines multiple (2+) fragments.
    """
    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent

        if node["type"] == "reaction" and depth <= 1:  # Check final and penultimate steps
            # This is the final or penultimate reaction
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Filter out any empty reactants
                valid_reactants = [r for r in reactants if r.strip()]

                # Check if there are 2+ significant reactants (more than 5 atoms each)
                significant_reactants = []
                for r in valid_reactants:
                    try:
                        mol = Chem.MolFromSmiles(r)
                        if mol and mol.GetNumAtoms() > 5:
                            significant_reactants.append(r)
                    except:
                        print(f"Warning: Could not parse reactant SMILES: {r}")

                # If there are 2 or more significant reactants, it's convergent
                if len(significant_reactants) >= 2:
                    # Additional check: ensure reactants are not just reagents
                    # by comparing their complexity to the product
                    product_mol = Chem.MolFromSmiles(product_part)
                    if product_mol:
                        product_atoms = product_mol.GetNumAtoms()
                        # Check if reactants contribute significantly to product
                        reactant_contribution = sum(
                            Chem.MolFromSmiles(r).GetNumAtoms()
                            for r in significant_reactants
                            if Chem.MolFromSmiles(r)
                        )
                        if (
                            reactant_contribution > 0.7 * product_atoms
                        ):  # At least 70% of product atoms come from reactants
                            is_convergent = True
                            print(
                                f"Found convergent synthesis with {len(significant_reactants)} significant fragments at depth {depth}: {significant_reactants}"
                            )
            else:
                print(f"Warning: Reaction node at depth {depth} missing metadata or rsmi")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Convergent synthesis detected: {is_convergent}")
    return is_convergent
