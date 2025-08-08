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
    This function detects a synthetic strategy involving Suzuki coupling
    to form biaryl bonds, typically in the middle stages of synthesis.
    """
    has_suzuki = False

    def dfs_traverse(node):
        nonlocal has_suzuki

        if node["type"] == "reaction":
            # Extract reaction SMILES
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            # Split into reactants and product
            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[2]

            # Check for Suzuki coupling patterns
            try:
                # Look for boronic acid and aryl bromide in reactants
                boronic_acid_pattern = Chem.MolFromSmarts("[c][B]([OH])[OH]")
                aryl_bromide_pattern = Chem.MolFromSmarts("[c][Br]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                has_boronic_acid = any(
                    mol and mol.HasSubstructMatch(boronic_acid_pattern) for mol in reactant_mols
                )
                has_aryl_bromide = any(
                    mol and mol.HasSubstructMatch(aryl_bromide_pattern) for mol in reactant_mols
                )

                # Check for biaryl in product
                product_mol = Chem.MolFromSmiles(product)
                biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

                if (
                    has_boronic_acid
                    and has_aryl_bromide
                    and product_mol
                    and product_mol.HasSubstructMatch(biaryl_pattern)
                ):
                    depth = node["metadata"].get("depth", -1)
                    print(f"Detected Suzuki coupling at depth {depth}")
                    has_suzuki = True
            except:
                pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_suzuki
