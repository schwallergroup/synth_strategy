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
    This function detects a linear synthesis strategy where an aromatic scaffold
    is preserved throughout the synthesis with sequential modifications.
    """
    # Track if we have a linear synthesis with preserved aromatic scaffold
    reaction_count = 0
    all_reactions_preserve_aromatic = True

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, all_reactions_preserve_aromatic

        if node["type"] == "reaction":
            reaction_count += 1

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if both reactant and product have aromatic rings
                main_reactant = reactants[0]  # Assuming first reactant is the main one
                reactant_mol = Chem.MolFromSmiles(main_reactant)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    aromatic_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1")

                    reactant_has_aromatic = reactant_mol.HasSubstructMatch(aromatic_pattern)
                    product_has_aromatic = product_mol.HasSubstructMatch(aromatic_pattern)

                    # If either doesn't have aromatic ring, the scaffold is not preserved
                    if not (reactant_has_aromatic and product_has_aromatic):
                        all_reactions_preserve_aromatic = False

                    # Check if this is a linear synthesis (one main reactant + small reagent)
                    if len(reactants) > 1:
                        # Check if second reactant is a small molecule (likely a reagent)
                        reagent_mol = Chem.MolFromSmiles(reactants[1])
                        if reagent_mol and reagent_mol.GetNumAtoms() > 6:  # If reagent is large
                            # This might be a convergent step rather than linear
                            print(f"Possible convergent step at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Strategy is present if we have multiple reactions that all preserve the aromatic scaffold
    result = reaction_count >= 2 and all_reactions_preserve_aromatic
    print(f"Linear synthesis with preserved aromatic scaffold detected: {result}")
    return result
