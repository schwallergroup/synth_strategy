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
    This function detects if the synthetic route involves multiple C-O bond formations,
    particularly focusing on aromatic C-O bond formations.
    """
    co_bond_formation_count = 0

    def dfs_traverse(node):
        nonlocal co_bond_formation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Convert SMILES to molecules
                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product = Chem.MolFromSmiles(product_smiles)

                    if all(r is not None for r in reactants) and product is not None:
                        # Check for phenol in reactants
                        phenol_pattern = Chem.MolFromSmarts("[OH][c]")
                        # Check for aryl ether in product
                        aryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[#6]")

                        reactants_have_phenol = any(
                            r.HasSubstructMatch(phenol_pattern) for r in reactants
                        )
                        product_has_aryl_ether = product.HasSubstructMatch(aryl_ether_pattern)

                        if reactants_have_phenol and product_has_aryl_ether:
                            print(f"C-O bond formation detected in reaction: {rsmi}")
                            co_bond_formation_count += 1
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if multiple C-O bond formations are detected
    return co_bond_formation_count >= 2
