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
    This function detects if the synthetic route involves coupling of
    thiophene and indole fragments.
    """
    # Track if thiophene-indole coupling is found
    coupling_found = False

    # Thiophene and indole SMARTS patterns
    thiophene_pattern = Chem.MolFromSmarts("c1cccs1")
    indole_pattern = Chem.MolFromSmarts("c1ccc2[n]ccc2c1")

    def dfs_traverse(node, depth=0):
        nonlocal coupling_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains both thiophene and indole
            product_mol = Chem.MolFromSmiles(product_smiles)
            if (
                product_mol
                and product_mol.HasSubstructMatch(thiophene_pattern)
                and product_mol.HasSubstructMatch(indole_pattern)
            ):

                # Check if reactants are separate thiophene and indole fragments
                thiophene_only = False
                indole_only = False

                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if not reactant_mol:
                        continue

                    has_thiophene = reactant_mol.HasSubstructMatch(thiophene_pattern)
                    has_indole = reactant_mol.HasSubstructMatch(indole_pattern)

                    if has_thiophene and not has_indole:
                        thiophene_only = True
                    elif has_indole and not has_thiophene:
                        indole_only = True

                if thiophene_only and indole_only:
                    coupling_found = True
                    print(f"Thiophene-indole coupling detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return coupling_found
