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
    Detects if the synthesis follows a linear (non-convergent) strategy.

    A linear synthesis has each reaction step using only one non-starting material
    (one molecule that needs to be synthesized). If any reaction step combines
    multiple complex molecules (non-starting materials), the synthesis is considered
    convergent.

    Args:
        route: A synthesis route tree following the SynthesisRoute schema

    Returns:
        bool: True if the synthesis is linear, False if convergent
    """
    is_linear = True
    debug = False  # Set to True to enable debug prints

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Extract reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = [r for r in rsmi.split(">")[0].split(".") if r.strip()]

            # Get non-starting material children that are used as reactants
            non_starting_material_children = []

            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", True):
                    # Check if this non-starting material is actually used in the reaction
                    # We need to compare the child SMILES with reactant SMILES
                    # This is challenging due to atom mapping differences, so we'll use a basic approach

                    # Get the unmapped SMILES for comparison
                    child_mol = Chem.MolFromSmiles(child["smiles"])
                    if child_mol:
                        child_smiles_unmapped = Chem.MolToSmiles(child_mol)

                        # Check if this child appears in any of the reactants
                        # We'll consider it a match if it's a substring of a reactant
                        # or if a reactant is a substring of it (accounting for atom mapping differences)
                        is_reactant = False
                        for reactant in reactants_smiles:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                reactant_smiles_unmapped = Chem.MolToSmiles(reactant_mol)

                                # Check if either is a substructure of the other
                                if (
                                    child_smiles_unmapped in reactant_smiles_unmapped
                                    or reactant_smiles_unmapped in child_smiles_unmapped
                                ):
                                    is_reactant = True
                                    break

                        if is_reactant:
                            non_starting_material_children.append(child)

            # If more than one non-starting material, it's convergent
            if len(non_starting_material_children) > 1:
                if debug:
                    print(
                        f"Depth {depth}: Convergent step detected with {len(non_starting_material_children)} non-starting materials"
                    )
                    print(f"Reaction SMILES: {rsmi}")
                    for i, child in enumerate(non_starting_material_children):
                        print(f"  Non-starting material {i+1}: {child['smiles']}")
                is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    return is_linear
