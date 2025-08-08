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
    Detects significant complexity increase through fragment coupling.
    Looks for reactions where multiple complex fragments are combined.
    """
    fragment_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal fragment_coupling_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we have multiple complex reactants
            complex_reactants = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.GetNumAtoms() > 10:  # Arbitrary threshold for "complex"
                    complex_reactants += 1

            # Check if product is significantly more complex than individual reactants
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and complex_reactants >= 2:
                # Check if product has more rings than any individual reactant
                product_rings = Chem.GetSSSR(product_mol)
                max_reactant_rings = 0

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        max_reactant_rings = max(max_reactant_rings, len(Chem.GetSSSR(mol)))

                if len(product_rings) > max_reactant_rings:
                    fragment_coupling_found = True
                    print(f"Found fragment coupling with complexity increase at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return fragment_coupling_found
