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
    This function detects a synthetic strategy involving sequential aromatic
    nucleophilic substitutions on a heterocyclic core.
    """
    substitution_count = 0

    def dfs_traverse(node):
        nonlocal substitution_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aromatic nucleophilic substitution (replacement of halogen with nitrogen)
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Look for aromatic carbon-nitrogen bonds in product
            c_n_pattern = Chem.MolFromSmarts("[c]-[#7]")

            # Look for aromatic carbon-halogen bonds in reactants
            c_x_pattern = Chem.MolFromSmarts("[c]-[#9,#17,#35,#53]")

            if product_mol and product_mol.HasSubstructMatch(c_n_pattern):
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                if any(mol and mol.HasSubstructMatch(c_x_pattern) for mol in reactant_mols):
                    substitution_count += 1
                    print(f"Aromatic nucleophilic substitution detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least one aromatic nucleophilic substitution
    if substitution_count >= 1:
        print(
            f"Strategy detected: Aromatic nucleophilic substitution sequence with {substitution_count} substitutions"
        )
        return True

    return False
