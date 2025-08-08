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
    This function detects a synthetic strategy involving nitro group reduction to amine.
    """
    nitro_to_amine_found = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant contains nitro group
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and any(r for r in reactant_mols if r):
                    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2][#6]")

                    # Check if product has amine but no nitro
                    product_has_amine = product_mol.HasSubstructMatch(amine_pattern)
                    product_has_nitro = product_mol.HasSubstructMatch(nitro_pattern)

                    # Check if any reactant has nitro but no amine
                    reactant_has_nitro = any(
                        r.HasSubstructMatch(nitro_pattern) for r in reactant_mols if r
                    )

                    if product_has_amine and not product_has_nitro and reactant_has_nitro:
                        print("Detected nitro reduction to amine")
                        nitro_to_amine_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return nitro_to_amine_found
