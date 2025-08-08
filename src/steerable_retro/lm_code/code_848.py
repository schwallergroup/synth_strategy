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
    This function detects nitro group reduction to amine in the synthetic route.
    """
    nitro_to_amine_detected = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create RDKit mol objects
            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # Check for nitro group in reactants
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                # Check for amine group in product
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                nitro_in_reactants = any(
                    mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols if mol
                )
                amine_in_product = product_mol and product_mol.HasSubstructMatch(amine_pattern)

                if nitro_in_reactants and amine_in_product:
                    print("Detected nitro to amine reduction")
                    nitro_to_amine_detected = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return nitro_to_amine_detected
