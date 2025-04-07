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
    Detects if the synthesis follows a linear strategy without convergent steps.

    A linear synthesis strategy is characterized by:
    1. Each reaction step has one main reactant that contributes most atoms to the product
    2. The synthesis follows a clear path from starting materials to the target molecule
    3. There are at least 2 reaction steps

    Returns:
        bool: True if the synthesis follows a linear strategy, False otherwise
    """
    # Common reagents that shouldn't count toward convergence
    common_reagents = [
        "CC(C)=O",  # acetone
        "O",  # water
        "O=[Mn](=O)(=O)O",  # permanganate
        "[K+]",  # potassium ion
        "C1COCCO1",  # dioxane
        "CC(C)(C)O",  # tert-butanol
        "[Na+]",  # sodium ion
        "[O-][I+3]([O-])([O-])O",  # periodate
        "CN1CCCC1=O",  # N-methylpyrrolidone
        "CC(=O)O",  # acetic acid
        "CCO",  # ethanol
        "CO",  # methanol
        "CC#N",  # acetonitrile
        "[H][H]",  # hydrogen
        "ClC(Cl)Cl",  # chloroform
        "C1CCCC1",  # cyclopentane (common solvent)
    ]

    # Initialize variables to track synthesis characteristics
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants and product from atom-mapped reaction SMILES
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Convert product to molecule
                product_mol = Chem.MolFromSmiles(product_part)
                if not product_mol:
                    print(f"Warning: Could not parse product SMILES: {product_part}")
                    return

                # Filter out common reagents and small molecules
                significant_reactants = []
                for r in reactants:
                    if not r:
                        continue

                    # Skip common reagents
                    is_common_reagent = False
                    for common in common_reagents:
                        if Chem.MolToSmiles(Chem.MolFromSmiles(r)) == Chem.MolToSmiles(
                            Chem.MolFromSmiles(common)
                        ):
                            is_common_reagent = True
                            break

                    if is_common_reagent:
                        continue

                    # Parse reactant molecule
                    reactant_mol = Chem.MolFromSmiles(r)
                    if not reactant_mol:
                        continue

                    # Skip very small molecules (likely reagents)
                    if reactant_mol.GetNumHeavyAtoms() <= 2:
                        continue

                    significant_reactants.append(r)

                # Check for convergent synthesis (multiple significant reactants)
                if len(significant_reactants) > 1:
                    # Count mapped atoms from each reactant to the product
                    atom_contributions = {}

                    for idx, reactant in enumerate(significant_reactants):
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if not reactant_mol:
                            continue

                        # Count mapped atoms
                        mapped_atom_count = 0
                        for atom in reactant_mol.GetAtoms():
                            if (
                                atom.GetProp("molAtomMapNumber")
                                if atom.HasProp("molAtomMapNumber")
                                else None
                            ):
                                mapped_atom_count += 1

                        atom_contributions[idx] = mapped_atom_count

                    # If multiple reactants contribute significantly, it's convergent
                    significant_contributors = 0
                    for idx, count in atom_contributions.items():
                        if count >= 3:  # At least 3 atoms contributed
                            significant_contributors += 1

                    if significant_contributors > 1:
                        print(f"Convergent step detected in reaction: {rsmi}")
                        print(f"Multiple significant reactants: {significant_reactants}")
                        is_linear = False

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # A linear synthesis should have at least 2 reactions
    strategy_present = is_linear and reaction_count >= 2

    if strategy_present:
        print(f"Linear synthesis strategy detected with {reaction_count} reactions")
    else:
        if not is_linear:
            print("Convergent steps detected, not a linear synthesis")
        else:
            print(f"Not enough reactions ({reaction_count}) to confirm linear synthesis strategy")

    return strategy_present
