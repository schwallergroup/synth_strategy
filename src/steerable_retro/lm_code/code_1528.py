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
    Detects if the synthesis route involves a C-C bond formation between
    two heterocyclic systems.
    """
    heterocycle_coupling_found = False

    def dfs_traverse(node):
        nonlocal heterocycle_coupling_found

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for heterocycles in reactants and a larger heterocyclic system in product
            heterocycle_pattern = Chem.MolFromSmarts("[a;!c]")  # Aromatic atom that's not carbon

            if heterocycle_pattern:
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]

                if product_mol and len(reactants_mols) >= 2:
                    # Count heterocycles in reactants and product
                    heterocycle_count_reactants = sum(
                        len(r.GetSubstructMatches(heterocycle_pattern)) for r in reactants_mols if r
                    )
                    heterocycle_count_product = (
                        len(product_mol.GetSubstructMatches(heterocycle_pattern))
                        if product_mol
                        else 0
                    )

                    # If product has heterocycles from both reactants, it's likely a coupling
                    if heterocycle_count_product > 0 and heterocycle_count_reactants >= 2:
                        heterocycle_coupling_found = True
                        print("Heterocycle coupling detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_coupling_found
