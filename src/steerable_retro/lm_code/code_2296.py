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
    Detects a synthetic strategy that utilizes hydrazide chemistry
    in multiple steps of the synthesis.
    """
    # Track hydrazide occurrences
    hydrazide_reactions = 0

    def dfs_traverse(node):
        nonlocal hydrazide_reactions

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not all(reactants) or not product:
                print("Warning: Could not parse all molecules in reaction")
                return

            # Check for hydrazide pattern in reactants or product
            hydrazide_pattern = Chem.MolFromSmarts("[NX3][NX3][CX3]=[OX1]")

            # Check product for hydrazide
            if product.HasSubstructMatch(hydrazide_pattern):
                hydrazide_reactions += 1
                print(f"Found hydrazide formation in product")

            # Check if hydrazine is a reactant
            hydrazine_pattern = Chem.MolFromSmarts("[NX3][NX3]")
            for reactant in reactants:
                if reactant.HasSubstructMatch(hydrazine_pattern) and not reactant.HasSubstructMatch(
                    hydrazide_pattern
                ):
                    hydrazide_reactions += 1
                    print(f"Found hydrazine as reactant")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple hydrazide-related reactions are found
    return hydrazide_reactions >= 2
