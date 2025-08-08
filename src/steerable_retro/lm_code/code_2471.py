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
    This function detects a synthetic strategy involving multiple ring transformations
    (at least two ring openings or formations).
    """
    ring_transformations = 0

    def dfs_traverse(node):
        nonlocal ring_transformations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                try:
                    # Process all reactants
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    reactant_mols = [mol for mol in reactant_mols if mol]  # Filter out None values

                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mols and product_mol:
                        # Count total rings in all reactants
                        total_reactant_rings = sum(len(Chem.GetSSSR(mol)) for mol in reactant_mols)
                        product_rings = len(Chem.GetSSSR(product_mol))

                        print(f"Reaction: {rsmi}")
                        print(
                            f"Reactant rings: {total_reactant_rings}, Product rings: {product_rings}"
                        )

                        # If ring count changes (either opening or formation), count it as a transformation
                        if product_rings != total_reactant_rings:
                            ring_transformations += 1
                            print(f"Found ring transformation, total: {ring_transformations}")
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if at least two ring transformations are detected
    return ring_transformations >= 2
