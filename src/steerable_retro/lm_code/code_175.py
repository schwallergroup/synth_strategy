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
    Detects if the synthesis route involves heterocycle formation via cyclization.
    Looks for an increase in the number of rings, specifically nitrogen-containing rings.
    """
    has_cyclization = False

    def dfs_traverse(node):
        nonlocal has_cyclization

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            product_mol = Chem.MolFromSmiles(product)
            if not product_mol:
                return

            # Count rings in product
            product_rings = Chem.GetSSSR(product_mol)
            product_n_rings = 0
            for ring in product_rings:
                for atom_idx in ring:
                    atom = product_mol.GetAtomWithIdx(atom_idx)
                    if atom.GetSymbol() == "N":
                        product_n_rings += 1
                        break

            # Count rings in reactants
            reactant_n_rings = 0
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                reactant_rings = Chem.GetSSSR(reactant_mol)
                for ring in reactant_rings:
                    for atom_idx in ring:
                        atom = reactant_mol.GetAtomWithIdx(atom_idx)
                        if atom.GetSymbol() == "N":
                            reactant_n_rings += 1
                            break

            # Check if nitrogen-containing rings increased
            if product_n_rings > reactant_n_rings:
                has_cyclization = True
                print(
                    f"Found heterocycle formation: {reactant_n_rings} → {product_n_rings} N-rings"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_cyclization
