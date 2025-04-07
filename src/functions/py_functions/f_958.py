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
    Detects if the synthetic route follows a linear synthesis approach
    (as opposed to convergent synthesis).

    Linear synthesis: Each reaction step adds one significant building block
    Convergent synthesis: At least one step combines multiple significant building blocks
    """
    is_linear = True

    def is_significant_reactant(smiles):
        """Helper function to determine if a reactant is significant (building block)
        rather than a reagent, catalyst, or solvent"""

        # Skip empty or invalid SMILES
        if not smiles or len(smiles) < 3:
            return False

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Count heavy atoms (non-hydrogen)
        heavy_atom_count = mol.GetNumHeavyAtoms()

        # Common reagents to exclude (regardless of size)
        common_reagents = [
            # Protecting groups and their reagents
            "CC(C)(C)OC(=O)",  # Boc anhydride or Boc-containing reagents
            "CC(=O)O",  # Acetate
            "CS(=O)(=O)",  # Mesyl group
            "CC(C)(C)Si",  # TMS or TBDMS groups
            # Common bases
            "CN(C)C",  # TEA, DIPEA, etc.
            "CC(C)N",  # Amines
            "CN",  # Simple amines
            # Common acids
            "OS(=O)(=O)O",  # Sulfuric acid
            "O=C(O)",  # Carboxylic acids
            # Reducing/oxidizing agents
            "B",  # Boron-containing (NaBH4, etc.)
            "OS(=O)",  # DMSO
            "O=N",  # Nitro compounds
            "Cl",  # Chlorides
            "Br",  # Bromides
            "I",  # Iodides
            # Coupling reagents
            "P(=O)",  # Phosphorus reagents
            "N=[N+]=[N-]",  # Azide
            # Solvents
            "CCO",  # Ethanol
            "CO",  # Methanol
            "CC(=O)",  # Acetone
            "ClCCl",  # DCM
            "c1ccccc1",  # Benzene
        ]

        # Check if the molecule contains any common reagent patterns
        for reagent in common_reagents:
            if reagent in smiles:
                return False

        # Consider molecules with more than 8 heavy atoms as significant
        # This threshold can be adjusted based on the specific chemistry
        if heavy_atom_count > 8:
            return True

        # For smaller molecules, check if they're likely building blocks
        # rather than reagents by looking at functional groups

        # If it has complex structure but small size, it might still be a building block
        if heavy_atom_count > 5 and (
            "c1" in smiles or "C1" in smiles
        ):  # Contains a ring
            return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if not is_linear:  # Early return if we already know it's not linear
            return

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count significant reactants (excluding reagents)
            significant_reactants = []
            for r in reactants_smiles:
                if is_significant_reactant(r):
                    significant_reactants.append(r)

            # If more than one significant reactant, it's not a linear synthesis
            if len(significant_reactants) > 1:
                is_linear = False
                print(
                    f"Found convergent step with multiple significant reactants: {significant_reactants}"
                )

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return is_linear
