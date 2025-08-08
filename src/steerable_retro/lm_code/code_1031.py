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
    This function detects if the synthesis involves multiple nitro-amine interconversions
    (nitro reduction or amine nitration).
    """
    interconversion_count = 0

    def dfs_traverse(node):
        nonlocal interconversion_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check for nitro reduction (NO2 to NH2)
                    if reactants_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c]-[N+](=[O])[O-]")
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[NH2]")):
                        interconversion_count += 1
                        print(f"Detected nitro reduction, count: {interconversion_count}")

                    # Check for amine nitration (NH2 to NO2)
                    if reactants_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c]-[NH2]")
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[N+](=[O])[O-]")):
                        interconversion_count += 1
                        print(f"Detected amine nitration, count: {interconversion_count}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Consider it a strategy if there are at least 2 interconversions
    multiple_interconversions = interconversion_count >= 2

    if multiple_interconversions:
        print(f"Multiple nitro-amine interconversions found: {interconversion_count}")

    return multiple_interconversions
